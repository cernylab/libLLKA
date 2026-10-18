"""Reading, editing and writing of mmCIF data (``llka_minicif.h``).

:func:`cif_file_to_structure` and :func:`cif_text_to_structure` import a
:class:`~llka.structure.Structure` from mmCIF. :class:`CifData` gives access to
the complete Cif data: blocks contain categories, categories contain items and
items contain values.

Adding items or values marks the data as *tainted*. Tainted data must be
validated with :meth:`CifData.detaint` before it can be written out.
"""

from __future__ import annotations

import dataclasses
import enum
import os
import weakref
from collections.abc import Iterable, Iterator
from ctypes import POINTER, byref, c_void_p
from typing import NamedTuple

from . import _capi
from ._capi import (
    LLKA_CifData,
    LLKA_CifDataBlock,
    LLKA_CifDataValue,
    LLKA_ImportedStructure,
    lib,
)
from .enums import RetCode
from .errors import BadDataError, CifParseError, InvalidArgumentError, error_for_retcode
from .structure import Structure

__all__ = [
    'CifData',
    'CifDataBlock',
    'CifDataCategory',
    'CifDataItem',
    'CifValue',
    'CifValueState',
    'ImportedStructure',
    'cif_data_to_string',
    'cif_file_to_data',
    'cif_file_to_structure',
    'cif_text_to_data',
    'cif_text_to_structure',
]

# Import options. LLKA_MINICIF_NORMALIZE is intentionally not exposed, it does not order atoms correctly.
_GET_CIFDATA = 1 << 1
_ALLOW_NO_ENTRY_CATEGORY = 1 << 2
_ALLOW_BROKEN_ATOMSITE = 1 << 3


class CifValueState(enum.IntEnum):
    """State of a Cif value (``LLKA_CifDataValueState``)."""

    SET = 0      # The value has a text
    NONE = 1     # Inapplicable value, written as "."
    UNKNOWN = 2  # Unknown value, written as "?"


class CifValue(NamedTuple):
    """A Cif value. ``text`` is ``None`` unless the state is :attr:`CifValueState.SET`."""

    text: str | None
    state: CifValueState = CifValueState.SET

    @classmethod
    def none(cls) -> CifValue:
        """Inapplicable value (``.``)."""
        return cls(None, CifValueState.NONE)

    @classmethod
    def unknown(cls) -> CifValue:
        """Unknown value (``?``)."""
        return cls(None, CifValueState.UNKNOWN)

    @classmethod
    def _from_c(cls, c: LLKA_CifDataValue) -> CifValue:
        state = CifValueState(c.state)
        return cls(_capi.decode(c.text) if state == CifValueState.SET else None, state)


CifValueLike = CifValue | str


def _value_to_c(value: CifValueLike, c: LLKA_CifDataValue) -> None:
    if isinstance(value, str):
        value = CifValue(value)
    elif not isinstance(value, CifValue):
        raise TypeError(f'Expected str or CifValue, got {type(value).__name__}')

    state = CifValueState(value.state)
    if state == CifValueState.SET:
        if value.text is None:
            raise ValueError('A Cif value in the SET state must have a text')
        c.text = _capi.encode(value.text)
    else:
        c.text = None
    c.state = state


def _values_to_c(values: Iterable[CifValueLike]):
    values = list(values)
    arr = (LLKA_CifDataValue * len(values))()
    for idx, v in enumerate(values):
        _value_to_c(v, arr[idx])
    return arr, len(values)


def _destroy_cif_data(ptr) -> None:
    lib.LLKA_destroyCifData(ptr)


class CifData:
    """Complete Cif data, i.e. a list of data blocks (``LLKA_CifData``)."""

    __slots__ = ('_ptr', '_generation', '_keepalive', '_finalizer', '__weakref__')

    def __init__(self) -> None:
        """Creates empty Cif data."""
        self._own(lib.LLKA_cifData_empty())

    def _own(self, ptr) -> None:
        self._ptr = ptr
        self._generation = 0
        self._keepalive = None
        self._finalizer = weakref.finalize(self, _destroy_cif_data, ptr)

    @classmethod
    def _take(cls, ptr) -> CifData:
        data = cls.__new__(cls)
        data._own(ptr)
        return data

    @classmethod
    def _borrow(cls, ptr, keepalive: object) -> CifData:
        data = cls.__new__(cls)
        data._ptr = ptr
        data._generation = 0
        data._keepalive = keepalive
        data._finalizer = None
        return data

    @classmethod
    def from_file(cls, path: str | os.PathLike[str]) -> CifData:
        """Same as :func:`cif_file_to_data`."""
        return cif_file_to_data(path)

    @classmethod
    def from_text(cls, text: str) -> CifData:
        """Same as :func:`cif_text_to_data`."""
        return cif_text_to_data(text)

    @property
    def _c(self) -> LLKA_CifData:
        return self._ptr.contents

    @property
    def blocks(self) -> list[CifDataBlock]:
        return [CifDataBlock(self, idx) for idx in range(self._c.nBlocks)]

    def __len__(self) -> int:
        return self._c.nBlocks

    def __getitem__(self, idx: int) -> CifDataBlock:
        n = self._c.nBlocks
        if idx < 0:
            idx += n
        if not 0 <= idx < n:
            raise IndexError('Block index out of range')
        return CifDataBlock(self, idx)

    def __iter__(self) -> Iterator[CifDataBlock]:
        return iter(self.blocks)

    def add_block(self, name: str = '') -> CifDataBlock:
        """Appends a new block. The name may be empty."""
        lib.LLKA_cifData_addBlock(self._ptr, _capi.encode(name))
        return CifDataBlock(self, self._c.nBlocks - 1)

    def delete_block(self, block: CifDataBlock | int) -> None:
        """Deletes a block. Wrappers of all blocks of this data become invalid."""
        idx = block._index if isinstance(block, CifDataBlock) else block
        if isinstance(block, CifDataBlock):
            block._check()
        _capi.check(lib.LLKA_cifData_deleteBlock(self._ptr, idx), f'Cannot delete block {idx}')
        self._generation += 1

    @property
    def is_tainted(self) -> bool:
        """Whether the data has been modified and needs to be validated with :meth:`detaint`."""
        return bool(lib.LLKA_cifData_isTainted(self._ptr))

    def detaint(self) -> None:
        """Validates the data and clears the taint flag.

        :raises BadDataError: A category has no items.
        :raises MismatchingSizesError: Items of a category do not have the same number of values.
        """
        if not self.is_tainted:
            return
        for block in self.blocks:
            for cat in block:
                if not cat._ptr().contents.firstItem:
                    raise BadDataError(f'Category {cat.name!r} has no items', RetCode.E_BAD_DATA)
        _capi.check(lib.LLKA_cifData_detaint(self._ptr), 'Cif data is not valid')

    def to_string(self, pretty: bool = False) -> str:
        """Same as :func:`cif_data_to_string`."""
        return cif_data_to_string(self, pretty)

    def __repr__(self) -> str:
        return f'<CifData with {len(self)} blocks>'


class CifDataBlock:
    """A data block of :class:`CifData`. Iterating over a block yields its categories."""

    __slots__ = ('_data', '_index', '_generation')

    def __init__(self, data: CifData, index: int) -> None:
        self._data = data
        self._index = index
        self._generation = data._generation

    def _check(self) -> None:
        if self._generation != self._data._generation:
            raise RuntimeError('Cif block is no longer valid because a block was deleted')

    def _ptr(self):
        self._check()
        return POINTER(LLKA_CifDataBlock)(self._data._c.blocks[self._index])

    @property
    def data(self) -> CifData:
        return self._data

    @property
    def name(self) -> str:
        return _capi.decode(self._ptr().contents.name)

    def _first(self):
        return self._ptr().contents.firstCategory

    def __iter__(self) -> Iterator[CifDataCategory]:
        cat = self._first()
        while cat:
            yield CifDataCategory(self, _capi.decode(cat.contents.name))
            cat = lib.LLKA_cifDataBlock_nextCategory(cat)

    @property
    def categories(self) -> list[CifDataCategory]:
        return list(self)

    def __len__(self) -> int:
        return sum(1 for _ in self)

    def find_category(self, name: str) -> CifDataCategory | None:
        """Returns the category with the given name or ``None``."""
        ptr = lib.LLKA_cifDataBlock_findCategory(self._ptr(), _capi.encode(name))
        return CifDataCategory(self, name) if ptr else None

    def __getitem__(self, name: str) -> CifDataCategory:
        cat = self.find_category(name)
        if cat is None:
            raise KeyError(name)
        return cat

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and self.find_category(name) is not None

    def add_category(self, name: str) -> CifDataCategory:
        """Appends a new, empty category.

        :raises InvalidArgumentError: The name is empty or the block already contains a category with the name.
        """
        ptr = lib.LLKA_cifDataBlock_addCategory(self._ptr(), _capi.encode(name))
        if not ptr:
            raise InvalidArgumentError(f'Cannot add category {name!r}: the name is empty or already used')
        return CifDataCategory(self, name)

    def delete_category(self, name: str) -> None:
        """Deletes a category.

        :raises InvalidArgumentError: There is no category with the name.
        """
        _capi.check(lib.LLKA_cifDataBlock_deleteCategory(self._ptr(), _capi.encode(name)), f'Cannot delete category {name!r}')

    def __delitem__(self, name: str) -> None:
        if name not in self:
            raise KeyError(name)
        self.delete_category(name)

    def __repr__(self) -> str:
        return f'<CifDataBlock {self.name!r}>'


class CifDataCategory:
    """A category of a Cif block. Iterating over a category yields its items."""

    __slots__ = ('_block', '_name')

    def __init__(self, block: CifDataBlock, name: str) -> None:
        self._block = block
        self._name = name

    def _ptr(self):
        ptr = lib.LLKA_cifDataBlock_findCategory(self._block._ptr(), _capi.encode(self._name))
        if not ptr:
            raise KeyError(f'Category {self._name!r} no longer exists')
        return ptr

    @property
    def block(self) -> CifDataBlock:
        return self._block

    @property
    def name(self) -> str:
        return self._name

    @property
    def is_loop(self) -> bool:
        """Whether the category is a loop. Updated by :meth:`CifData.detaint`."""
        return bool(lib.LLKA_cifDataCategory_isLoop(self._ptr()))

    @property
    def next(self) -> CifDataCategory | None:
        """The following category in the block."""
        ptr = lib.LLKA_cifDataBlock_nextCategory(self._ptr())
        return CifDataCategory(self._block, _capi.decode(ptr.contents.name)) if ptr else None

    @property
    def previous(self) -> CifDataCategory | None:
        """The preceding category in the block."""
        ptr = lib.LLKA_cifDataBlock_previousCategory(self._ptr())
        return CifDataCategory(self._block, _capi.decode(ptr.contents.name)) if ptr else None

    def __iter__(self) -> Iterator[CifDataItem]:
        item = self._ptr().contents.firstItem
        while item:
            yield CifDataItem(self, _capi.decode(item.contents.keyword))
            item = lib.LLKA_cifDataCategory_nextItem(item)

    @property
    def items(self) -> list[CifDataItem]:
        return list(self)

    def __len__(self) -> int:
        return sum(1 for _ in self)

    def find_item(self, keyword: str) -> CifDataItem | None:
        """Returns the item with the given keyword or ``None``."""
        ptr = lib.LLKA_cifDataCategory_findItem(self._ptr(), _capi.encode(keyword))
        return CifDataItem(self, keyword) if ptr else None

    def __getitem__(self, keyword: str) -> CifDataItem:
        item = self.find_item(keyword)
        if item is None:
            raise KeyError(keyword)
        return item

    def __contains__(self, keyword: object) -> bool:
        return isinstance(keyword, str) and self.find_item(keyword) is not None

    def add_item(self, keyword: str, values: Iterable[CifValueLike] | None = None) -> CifDataItem:
        """Appends a new item, optionally with values. Taints the Cif data.

        :raises InvalidArgumentError: The keyword is empty or the category already contains an item with the keyword.
        """
        ptr = lib.LLKA_cifDataCategory_addItem(self._ptr(), _capi.encode(keyword))
        if not ptr:
            raise InvalidArgumentError(f'Cannot add item {keyword!r}: the keyword is empty or already used')
        item = CifDataItem(self, keyword)
        if values is not None:
            item.set_values(values)
        return item

    def delete_item(self, keyword: str) -> None:
        """Deletes an item.

        :raises InvalidArgumentError: There is no item with the keyword.
        """
        _capi.check(lib.LLKA_cifDataCategory_deleteItem(self._ptr(), _capi.encode(keyword)), f'Cannot delete item {keyword!r}')

    def __delitem__(self, keyword: str) -> None:
        if keyword not in self:
            raise KeyError(keyword)
        self.delete_item(keyword)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, CifDataCategory):
            return self._block._data is other._block._data and self._block._index == other._block._index and self._name == other._name
        return NotImplemented

    def __hash__(self) -> int:
        return hash((id(self._block._data), self._block._index, self._name))

    def __repr__(self) -> str:
        return f'<CifDataCategory {self._name!r}>'


class CifDataItem:
    """An item of a Cif category, i.e. a keyword and a list of values."""

    __slots__ = ('_category', '_keyword')

    def __init__(self, category: CifDataCategory, keyword: str) -> None:
        self._category = category
        self._keyword = keyword

    def _ptr(self):
        ptr = lib.LLKA_cifDataCategory_findItem(self._category._ptr(), _capi.encode(self._keyword))
        if not ptr:
            raise KeyError(f'Item {self._keyword!r} no longer exists')
        return ptr

    @property
    def category(self) -> CifDataCategory:
        return self._category

    @property
    def keyword(self) -> str:
        return self._keyword

    @property
    def next(self) -> CifDataItem | None:
        """The following item in the category."""
        ptr = lib.LLKA_cifDataCategory_nextItem(self._ptr())
        return CifDataItem(self._category, _capi.decode(ptr.contents.keyword)) if ptr else None

    @property
    def previous(self) -> CifDataItem | None:
        """The preceding item in the category."""
        ptr = lib.LLKA_cifDataCategory_previousItem(self._ptr())
        return CifDataItem(self._category, _capi.decode(ptr.contents.keyword)) if ptr else None

    @property
    def values(self) -> list[CifValue]:
        c = self._ptr().contents
        return [CifValue._from_c(c.values[idx]) for idx in range(c.nValues)]

    @values.setter
    def values(self, values: Iterable[CifValueLike]) -> None:
        self.set_values(values)

    def __len__(self) -> int:
        return self._ptr().contents.nValues

    def __getitem__(self, idx: int) -> CifValue:
        c = self._ptr().contents
        if idx < 0:
            idx += c.nValues
        if not 0 <= idx < c.nValues:
            raise IndexError('Cif value index out of range')
        return CifValue._from_c(c.values[idx])

    def __iter__(self) -> Iterator[CifValue]:
        return iter(self.values)

    def add_value(self, value: CifValueLike) -> None:
        """Appends a value. Taints the Cif data."""
        c = LLKA_CifDataValue()
        _value_to_c(value, c)
        lib.LLKA_cifDataItem_addValue(self._ptr(), byref(c))

    def add_values(self, values: Iterable[CifValueLike]) -> None:
        """Appends values. Taints the Cif data."""
        arr, n = _values_to_c(values)
        lib.LLKA_cifDataItem_addValues(self._ptr(), arr, n)

    def set_values(self, values: Iterable[CifValueLike]) -> None:
        """Replaces all values. Taints the Cif data."""
        arr, n = _values_to_c(values)
        lib.LLKA_cifDataItem_setValues(self._ptr(), arr, n)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, CifDataItem):
            return self._category == other._category and self._keyword == other._keyword
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self._category, self._keyword))

    def __repr__(self) -> str:
        return f'<CifDataItem {self._keyword!r}>'


@dataclasses.dataclass(frozen=True, slots=True)
class ImportedStructure:
    """Structure imported from mmCIF.

    ``cif_data`` contains the complete Cif data if it was requested, otherwise ``None``.
    """

    entry_id: str
    structure: Structure
    cif_data: CifData | None


class _ImportedStructureOwner:
    __slots__ = ('_c', '_finalizer', '__weakref__')

    def __init__(self, c: LLKA_ImportedStructure) -> None:
        self._c = c
        self._finalizer = weakref.finalize(self, lib.LLKA_destroyImportedStructure, byref(c))


def _options(get_cif_data: bool, allow_no_entry_category: bool, allow_broken_atom_site: bool) -> int:
    options = 0
    if get_cif_data:
        options |= _GET_CIFDATA
    if allow_no_entry_category:
        options |= _ALLOW_NO_ENTRY_CATEGORY
    if allow_broken_atom_site:
        options |= _ALLOW_BROKEN_ATOMSITE
    return options


def _check_parse(code: int, error: c_void_p, context: str) -> None:
    message = _capi.take_string(error)
    if code == RetCode.OK:
        return
    if code == RetCode.E_BAD_DATA and message:
        raise CifParseError(message, RetCode.E_BAD_DATA)
    raise error_for_retcode(code, f'{context}: {_capi.error_to_string(code)}' + (f' ({message})' if message else ''))


def _imported_structure(c: LLKA_ImportedStructure) -> ImportedStructure:
    owner = _ImportedStructureOwner(c)
    cif_data = CifData._borrow(c.cifData, owner) if c.cifData else None
    return ImportedStructure(
        entry_id=_capi.decode(c.entry.id),
        structure=Structure._borrow(c.structure, owner),
        cif_data=cif_data,
    )


def cif_file_to_structure(
    path: str | os.PathLike[str],
    *,
    get_cif_data: bool = False,
    allow_no_entry_category: bool = False,
    allow_broken_atom_site: bool = False,
) -> ImportedStructure:
    """Imports a structure from an mmCIF file.

    :param get_cif_data: Also return the complete Cif data.
    :param allow_no_entry_category: Tolerate Cif data without the ``entry`` category.
    :param allow_broken_atom_site: Tolerate ``atom_site`` entries with missing ``label_`` fields
                                   and fill them in from the ``auth_`` fields.
    :raises NoFileError: The file does not exist.
    :raises NoDataError: The file is empty.
    :raises CannotReadFileError: The file cannot be read.
    :raises CifParseError: The file is malformed.
    """
    fn = _capi.require_function('LLKA_cifFileToStructure')
    c = LLKA_ImportedStructure()
    error = c_void_p()
    code = fn(
        _capi.path_arg(path), byref(c), byref(error),
        _options(get_cif_data, allow_no_entry_category, allow_broken_atom_site),
    )
    _check_parse(code, error, f'Cannot import structure from {os.fspath(path)}')
    return _imported_structure(c)


def cif_text_to_structure(
    text: str,
    *,
    get_cif_data: bool = False,
    allow_no_entry_category: bool = False,
    allow_broken_atom_site: bool = False,
) -> ImportedStructure:
    """Imports a structure from mmCIF text. See :func:`cif_file_to_structure` for the options.

    :raises CifParseError: The text is malformed.
    """
    c = LLKA_ImportedStructure()
    error = c_void_p()
    code = lib.LLKA_cifTextToStructure(
        _capi.encode(text), byref(c), byref(error),
        _options(get_cif_data, allow_no_entry_category, allow_broken_atom_site),
    )
    _check_parse(code, error, 'Cannot import structure')
    return _imported_structure(c)


def cif_file_to_data(path: str | os.PathLike[str]) -> CifData:
    """Reads Cif data from a file.

    :raises NoFileError: The file does not exist.
    :raises NoDataError: The file is empty.
    :raises CannotReadFileError: The file cannot be read.
    :raises CifParseError: The file is malformed.
    """
    fn = _capi.require_function('LLKA_cifFileToData')
    ptr = POINTER(LLKA_CifData)()
    error = c_void_p()
    code = fn(_capi.path_arg(path), byref(ptr), byref(error))
    _check_parse(code, error, f'Cannot read Cif data from {os.fspath(path)}')
    return CifData._take(ptr)


def cif_text_to_data(text: str) -> CifData:
    """Parses Cif data from text.

    :raises CifParseError: The text is malformed.
    """
    ptr = POINTER(LLKA_CifData)()
    error = c_void_p()
    code = lib.LLKA_cifTextToData(_capi.encode(text), byref(ptr), byref(error))
    _check_parse(code, error, 'Cannot parse Cif data')
    return CifData._take(ptr)


def cif_data_to_string(data: CifData, pretty: bool = False) -> str:
    """Writes Cif data as mmCIF text.

    :param pretty: Pad values to make the output easier to read for humans. Pretty output is larger.
    :raises BadDataError: The data is tainted or invalid.
    """
    out = c_void_p()
    _capi.check(lib.LLKA_cifDataToString(data._ptr, bool(pretty), byref(out)), 'Cannot write Cif data')
    return _capi.take_string(out) or ''

