"""Atoms and structures (``llka_structure.h``).

:class:`Atom` is an immutable value. :class:`Structure` is a mutable sequence of
atoms whose data lives in memory managed by libLLKA; atoms are converted to
Python objects only when they are accessed. :class:`StructureView` is a
read-only view of some atoms of a structure.
"""

from __future__ import annotations

import ctypes
import dataclasses
import weakref
from collections.abc import Iterable, Iterator, MutableSequence, Sequence
from ctypes import POINTER, byref
from typing import Any, overload

from . import _capi
from ._capi import LLKA_Atom, LLKA_Structure, LLKA_Structures, LLKA_StructureView, lib
from .core import Point, _to_point
from .enums import BaseKind
from .errors import InvalidArgumentError

__all__ = [
    'Atom',
    'Structure',
    'StructureLike',
    'StructureView',
    'base_kind',
    'compare_atoms',
]

_ATOM_PTR = POINTER(LLKA_Atom)
_ATOM_SIZE = ctypes.sizeof(LLKA_Atom)


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class Atom:
    """An atom. The fields correspond to the mmCIF ``atom_site`` category.

    ``auth_atom_id``, ``auth_comp_id`` and ``auth_asym_id`` default to the
    corresponding ``label_`` values. An empty ``label_alt_id`` means that the
    atom has no alternate position, an empty ``pdbx_PDB_ins_code`` means that
    there is no insertion code.
    """

    id: int
    type_symbol: str
    label_atom_id: str
    label_entity_id: str
    label_comp_id: str
    label_asym_id: str
    label_seq_id: int
    auth_seq_id: int
    coords: Point
    auth_atom_id: str = None  # type: ignore[assignment]
    auth_comp_id: str = None  # type: ignore[assignment]
    auth_asym_id: str = None  # type: ignore[assignment]
    pdbx_PDB_model_num: int = 1
    pdbx_PDB_ins_code: str = ''
    label_alt_id: str = ''

    def __post_init__(self) -> None:
        if self.auth_atom_id is None:
            object.__setattr__(self, 'auth_atom_id', self.label_atom_id)
        if self.auth_comp_id is None:
            object.__setattr__(self, 'auth_comp_id', self.label_comp_id)
        if self.auth_asym_id is None:
            object.__setattr__(self, 'auth_asym_id', self.label_asym_id)
        if not isinstance(self.coords, Point):
            object.__setattr__(self, 'coords', _to_point(self.coords))
        if len(self.label_alt_id) > 1 or not self.label_alt_id.isascii():
            raise ValueError('label_alt_id must be an empty string or a single ASCII character')

    def replace(self, **changes: Any) -> Atom:
        """Returns a copy of the atom with some fields replaced."""
        return dataclasses.replace(self, **changes)

    @classmethod
    def _from_c(cls, c: LLKA_Atom) -> Atom:
        decode = _capi.decode
        alt_id = c.label_alt_id
        coords = c.coords
        return cls(
            id=c.id,
            type_symbol=decode(c.type_symbol),
            label_atom_id=decode(c.label_atom_id),
            label_entity_id=decode(c.label_entity_id),
            label_comp_id=decode(c.label_comp_id),
            label_asym_id=decode(c.label_asym_id),
            label_seq_id=c.label_seq_id,
            auth_seq_id=c.auth_seq_id,
            coords=Point(coords.x, coords.y, coords.z),
            auth_atom_id=decode(c.auth_atom_id),
            auth_comp_id=decode(c.auth_comp_id),
            auth_asym_id=decode(c.auth_asym_id),
            pdbx_PDB_model_num=c.pdbx_PDB_model_num,
            pdbx_PDB_ins_code=decode(c.pdbx_PDB_ins_code),
            label_alt_id='' if alt_id == _capi.LLKA_NO_ALTID else alt_id.decode('ascii'),
        )

    def _fill_c(self, c: LLKA_Atom) -> None:
        """Fills a ctypes atom. The strings are owned by Python so the atom must not be passed to libLLKA destructors."""
        encode = _capi.encode
        c.type_symbol = encode(self.type_symbol)
        c.label_atom_id = encode(self.label_atom_id)
        c.label_entity_id = encode(self.label_entity_id)
        c.label_comp_id = encode(self.label_comp_id)
        c.label_asym_id = encode(self.label_asym_id)
        c.auth_atom_id = encode(self.auth_atom_id)
        c.auth_comp_id = encode(self.auth_comp_id)
        c.auth_asym_id = encode(self.auth_asym_id)
        c.coords.x, c.coords.y, c.coords.z = self.coords
        c.id = self.id
        c.label_seq_id = self.label_seq_id
        c.auth_seq_id = self.auth_seq_id
        c.pdbx_PDB_model_num = self.pdbx_PDB_model_num
        c.pdbx_PDB_ins_code = encode(self.pdbx_PDB_ins_code)
        c.label_alt_id = _alt_id_to_c(self.label_alt_id)

    def _to_c(self) -> LLKA_Atom:
        c = LLKA_Atom()
        self._fill_c(c)
        return c


def _alt_id_to_c(alt_id: str) -> bytes:
    return alt_id.encode('ascii') if alt_id else _capi.LLKA_NO_ALTID


def _atoms_to_c_array(atoms: Sequence[Atom]) -> ctypes.Array:
    arr = (LLKA_Atom * len(atoms))()
    for idx, atom in enumerate(atoms):
        if not isinstance(atom, Atom):
            raise TypeError(f'Expected Atom, got {type(atom).__name__}')
        atom._fill_c(arr[idx])
    return arr


def _destroy_structure(c: LLKA_Structure) -> None:
    lib.LLKA_destroyStructure(byref(c))


def _destroy_structures(c: LLKA_Structures) -> None:
    lib.LLKA_destroyStructures(byref(c))


class Structure(MutableSequence[Atom]):
    """A structure, i.e. a list of atoms stored in memory managed by libLLKA (``LLKA_Structure``).

    Indexing returns :class:`Atom` values, so modifying an atom requires
    assigning a new atom, e.g. ``stru[0] = stru[0].replace(label_seq_id=2)``.
    """

    __slots__ = ('_c', '_keepalive', '_generation', '_finalizer', '__weakref__')

    def __init__(self, atoms: Iterable[Atom] = ()) -> None:
        atoms = list(atoms)
        c = LLKA_Structure()
        lib.LLKA_initStructure(_atoms_to_c_array(atoms), len(atoms), byref(c))
        self._own(c)

    def _own(self, c: LLKA_Structure) -> None:
        self._c = c
        self._keepalive = None
        self._generation = 0
        self._finalizer = weakref.finalize(self, _destroy_structure, c)

    @classmethod
    def _take(cls, c: LLKA_Structure) -> Structure:
        """Takes ownership of a structure created by libLLKA."""
        stru = cls.__new__(cls)
        stru._own(c)
        return stru

    @classmethod
    def _borrow(cls, c: LLKA_Structure, keepalive: object) -> Structure:
        """Wraps a structure that is owned by another object, e.g. an element of an ``LLKA_Structures`` array."""
        stru = cls.__new__(cls)
        stru._c = c
        stru._keepalive = keepalive
        stru._generation = 0
        stru._finalizer = None
        return stru

    @classmethod
    def _take_array(cls, c: LLKA_Structures) -> list[Structure]:
        """Takes ownership of an array of structures created by libLLKA."""
        owner = _StructuresOwner(c)
        return [cls._borrow(c.strus[idx], owner) for idx in range(c.nStrus)]

    # Sequence protocol

    def __len__(self) -> int:
        return self._c.nAtoms

    def _index(self, idx: int) -> int:
        n = self._c.nAtoms
        if idx < 0:
            idx += n
        if not 0 <= idx < n:
            raise IndexError('Structure index out of range')
        return idx

    def _atom_ptr(self, idx: int):
        return ctypes.cast(ctypes.addressof(self._c.atoms.contents) + idx * _ATOM_SIZE, _ATOM_PTR)

    @overload
    def __getitem__(self, idx: int) -> Atom: ...
    @overload
    def __getitem__(self, idx: slice) -> Structure: ...

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            indices = range(*idx.indices(len(self)))
            ptrs = (_ATOM_PTR * len(indices))(*[self._atom_ptr(i) for i in indices])
            return Structure._take(lib.LLKA_makeStructureFromPtrs(ptrs, len(indices)))
        return Atom._from_c(self._c.atoms[self._index(idx)])

    def __iter__(self) -> Iterator[Atom]:
        c = self._c
        for idx in range(c.nAtoms):
            yield Atom._from_c(c.atoms[idx])

    def __setitem__(self, idx, value) -> None:
        if isinstance(idx, slice):
            indices = range(*idx.indices(len(self)))
            values = list(value)
            if idx.step not in (None, 1) or len(values) != len(indices):
                # Size-changing slice assignment, rebuild the whole structure
                atoms = list(self)
                atoms[idx] = values
                self._replace_with(atoms)
            else:
                for i, v in zip(indices, values):
                    self[i] = v
            return

        if not isinstance(value, Atom):
            raise TypeError(f'Expected Atom, got {type(value).__name__}')
        slot = self._c.atoms[self._index(idx)]
        lib.LLKA_destroyAtom(byref(slot))
        _init_c_atom(value, slot)

    def __delitem__(self, idx) -> None:
        n = len(self)
        if isinstance(idx, slice):
            removed = set(range(*idx.indices(n)))
        else:
            removed = {self._index(idx)}
        self._rebuild([self._atom_ptr(i) for i in range(n) if i not in removed])

    def insert(self, index: int, value: Atom) -> None:
        if not isinstance(value, Atom):
            raise TypeError(f'Expected Atom, got {type(value).__name__}')
        n = len(self)
        if index < 0:
            index = max(0, index + n)
        index = min(index, n)

        tmp = value._to_c()
        ptrs = [self._atom_ptr(i) for i in range(n)]
        ptrs.insert(index, ctypes.pointer(tmp))
        self._rebuild(ptrs)

    def _rebuild(self, ptrs: list) -> None:
        """Replaces the contents of the structure with deep copies of atoms pointed to by ``ptrs``."""
        arr = (_ATOM_PTR * len(ptrs))(*ptrs)
        new = lib.LLKA_makeStructureFromPtrs(arr, len(ptrs))
        self._swap_contents(new)

    def _replace_with(self, atoms: Sequence[Atom]) -> None:
        new = lib.LLKA_makeStructure(_atoms_to_c_array(atoms), len(atoms))
        self._swap_contents(new)

    def _swap_contents(self, new: LLKA_Structure) -> None:
        old = LLKA_Structure(self._c.atoms, self._c.nAtoms)
        self._c.atoms = new.atoms
        self._c.nAtoms = new.nAtoms
        _destroy_structure(old)
        self._invalidate_views()

    def _invalidate_views(self) -> None:
        self._generation += 1

    # Mutation

    def append(self, atom: Atom) -> None:
        """Appends a copy of an atom to the structure."""
        if not isinstance(atom, Atom):
            raise TypeError(f'Expected Atom, got {type(atom).__name__}')
        tmp = atom._to_c()
        lib.LLKA_appendAtom(byref(tmp), byref(self._c))
        self._invalidate_views()

    def extend(self, atoms: Iterable[Atom]) -> None:
        """Appends copies of atoms to the structure."""
        if isinstance(atoms, Structure):
            if atoms is self:
                atoms = atoms.copy()
            lib.LLKA_appendStructure(byref(atoms._c), byref(self._c))
        else:
            tmp = Structure(atoms)
            lib.LLKA_appendStructure(byref(tmp._c), byref(self._c))
        self._invalidate_views()

    def remove_atom_by_id(self, atom_id: int) -> None:
        """Removes the first atom with the given ``id``.

        :raises ValueError: There is no atom with the given ``id``.
        """
        n = len(self)
        lib.LLKA_removeAtomById(atom_id, byref(self._c))
        if len(self) == n:
            raise ValueError(f'Structure contains no atom with id {atom_id}')
        self._invalidate_views()

    def copy(self) -> Structure:
        """Returns a deep copy of the structure."""
        c = LLKA_Structure()
        lib.LLKA_duplicateStructure(byref(self._c), byref(c))
        return Structure._take(c)

    __copy__ = copy

    def __deepcopy__(self, memo: dict) -> Structure:
        return self.copy()

    # Queries

    def find_atom(
        self,
        label_atom_id: str | None = None,
        label_comp_id: str | None = None,
        label_asym_id: str | None = None,
        label_seq_id: int | None = None,
        label_alt_id: str = '',
        pdbx_PDB_ins_code: str = '',
        pdbx_PDB_model_num: int = 1,
    ) -> Atom | None:
        """Returns the first atom that matches the given criteria, or ``None``.

        Criteria that are ``None`` (or empty strings for ``label_alt_id`` and
        ``pdbx_PDB_ins_code``) are disregarded.
        """
        idx = self.find_atom_index(
            label_atom_id, label_comp_id, label_asym_id, label_seq_id, label_alt_id, pdbx_PDB_ins_code, pdbx_PDB_model_num
        )
        return None if idx is None else self[idx]

    def find_atom_index(
        self,
        label_atom_id: str | None = None,
        label_comp_id: str | None = None,
        label_asym_id: str | None = None,
        label_seq_id: int | None = None,
        label_alt_id: str = '',
        pdbx_PDB_ins_code: str = '',
        pdbx_PDB_model_num: int = 1,
    ) -> int | None:
        """Like :meth:`find_atom` but returns the index of the atom in the structure."""
        if len(self) == 0:
            return None

        ptr = lib.LLKA_findAtom(
            byref(self._c),
            *_find_atom_criteria(
                label_atom_id, label_comp_id, label_asym_id, label_seq_id, label_alt_id, pdbx_PDB_ins_code, pdbx_PDB_model_num
            ),
        )
        if not ptr:
            return None
        offset = ctypes.cast(ptr, ctypes.c_void_p).value - ctypes.addressof(self._c.atoms.contents)
        return offset // _ATOM_SIZE

    def split_by_alt_ids(self) -> list[tuple[str, Structure]]:
        """Splits the structure into substructures where each contains atoms with only one alternate position.

        Atoms without an alternate position are included in all substructures.

        :returns: List of ``(alt_id, substructure)`` pairs. If the structure has no alternate
                  positions, the list contains a single pair with an empty ``alt_id``.
        """
        alts = _capi.LLKA_AlternatePositions()
        c_strus = lib.LLKA_splitByAltIds(byref(self._c), byref(alts))
        try:
            alt_ids = [alts.positions[idx].decode('ascii') for idx in range(alts.nPositions)]
        finally:
            if alts.nPositions > 0:
                lib.LLKA_destroyAlternatePositions(byref(alts))
        strus = Structure._take_array(c_strus)
        if not alt_ids:
            alt_ids = [''] * len(strus)
        return list(zip(alt_ids, strus))

    def split_to_dinucleotide_steps(self) -> list[Structure]:
        """Splits the structure into dinucleotide steps.

        The structure must be ordered by model, chain and residue.

        :raises BadDataError: The structure cannot be split into steps.
        """
        c_strus = LLKA_Structures()
        _capi.check(
            lib.LLKA_splitStructureToDinucleotideSteps(byref(self._c), byref(c_strus)),
            'Cannot split structure to dinucleotide steps',
        )
        return Structure._take_array(c_strus)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Structure | StructureView):
            return len(self) == len(other) and all(a == b for a, b in zip(self, other))
        return NotImplemented

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f'<Structure with {len(self)} atoms>'


class _StructuresOwner:
    """Owns an ``LLKA_Structures`` array whose elements are shared by several :class:`Structure` objects."""

    __slots__ = ('_c', '_finalizer', '__weakref__')

    def __init__(self, c: LLKA_Structures) -> None:
        self._c = c
        self._finalizer = weakref.finalize(self, _destroy_structures, c)


def _find_atom_criteria(
    label_atom_id: str | None,
    label_comp_id: str | None,
    label_asym_id: str | None,
    label_seq_id: int | None,
    label_alt_id: str,
    pdbx_PDB_ins_code: str,
    pdbx_PDB_model_num: int,
) -> tuple:
    """Converts search criteria to the arguments of ``LLKA_findAtom()`` and ``LLKA_findAtomView()``."""
    def opt(s: str | None) -> bytes | None:
        return None if s is None else _capi.encode(s)

    return (
        opt(label_atom_id),
        opt(label_comp_id),
        opt(label_asym_id),
        -1 if label_seq_id is None else label_seq_id,
        _alt_id_to_c(label_alt_id),
        _capi.encode(pdbx_PDB_ins_code),
        pdbx_PDB_model_num,
    )


def _init_c_atom(atom: Atom, slot: LLKA_Atom) -> None:
    """Initializes an atom in memory owned by libLLKA with copies of the atom's data."""
    encode = _capi.encode
    coords = atom.coords._to_c()
    lib.LLKA_initAtom(
        atom.id,
        encode(atom.type_symbol),
        encode(atom.label_atom_id),
        encode(atom.label_entity_id),
        encode(atom.label_comp_id),
        encode(atom.label_asym_id),
        encode(atom.auth_atom_id),
        encode(atom.auth_comp_id),
        encode(atom.auth_asym_id),
        atom.label_seq_id,
        _alt_id_to_c(atom.label_alt_id),
        atom.auth_seq_id,
        encode(atom.pdbx_PDB_ins_code),
        atom.pdbx_PDB_model_num,
        byref(coords),
        byref(slot),
    )


class StructureView(Sequence[Atom]):
    """Read-only view of atoms of a :class:`Structure` (``LLKA_StructureView``).

    The view refers to the atoms of its parent structure. It becomes invalid
    when atoms are added to or removed from the parent structure; accessing
    an invalid view raises :class:`RuntimeError`.
    """

    __slots__ = ('_c', '_parent', '_generation', '_keepalive', '_finalizer', '__weakref__')

    @classmethod
    def _take(cls, c: LLKA_StructureView, parent: Structure) -> StructureView:
        """Takes ownership of a view created by libLLKA."""
        view = cls.__new__(cls)
        view._c = c
        view._parent = parent
        view._generation = parent._generation
        view._keepalive = None
        view._finalizer = weakref.finalize(view, lib.LLKA_destroyStructureView, byref(c))
        return view

    @classmethod
    def _borrow(cls, atoms, n_atoms: int, parent: Structure, keepalive: object) -> StructureView:
        """Wraps an array of atom pointers owned by another object."""
        view = cls.__new__(cls)
        view._c = LLKA_StructureView(atoms, n_atoms, n_atoms)
        view._parent = parent
        view._generation = parent._generation
        view._keepalive = keepalive
        view._finalizer = None
        return view

    @property
    def structure(self) -> Structure:
        """The structure that this view refers to."""
        return self._parent

    @property
    def is_valid(self) -> bool:
        """Whether the parent structure has not been modified in a way that invalidates the view."""
        return self._generation == self._parent._generation

    def _check(self) -> None:
        if not self.is_valid:
            raise RuntimeError('Structure view is no longer valid because its parent structure has been modified')

    def __len__(self) -> int:
        return self._c.nAtoms

    @overload
    def __getitem__(self, idx: int) -> Atom: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Atom]: ...

    def __getitem__(self, idx):
        self._check()
        if isinstance(idx, slice):
            return [self[i] for i in range(*idx.indices(len(self)))]
        n = self._c.nAtoms
        if idx < 0:
            idx += n
        if not 0 <= idx < n:
            raise IndexError('StructureView index out of range')
        return Atom._from_c(self._c.atoms[idx].contents)

    def __iter__(self) -> Iterator[Atom]:
        self._check()
        c = self._c
        for idx in range(c.nAtoms):
            yield Atom._from_c(c.atoms[idx].contents)

    def find_atom(
        self,
        label_atom_id: str | None = None,
        label_comp_id: str | None = None,
        label_asym_id: str | None = None,
        label_seq_id: int | None = None,
        label_alt_id: str = '',
        pdbx_PDB_ins_code: str = '',
        pdbx_PDB_model_num: int = 1,
    ) -> Atom | None:
        """Returns the first viewed atom that matches the given criteria, or ``None``.

        The criteria are the same as for :meth:`Structure.find_atom`.
        """
        idx = self.find_atom_index(
            label_atom_id, label_comp_id, label_asym_id, label_seq_id, label_alt_id, pdbx_PDB_ins_code, pdbx_PDB_model_num
        )
        return None if idx is None else self[idx]

    def find_atom_index(
        self,
        label_atom_id: str | None = None,
        label_comp_id: str | None = None,
        label_asym_id: str | None = None,
        label_seq_id: int | None = None,
        label_alt_id: str = '',
        pdbx_PDB_ins_code: str = '',
        pdbx_PDB_model_num: int = 1,
    ) -> int | None:
        """Like :meth:`find_atom` but returns the index of the atom in the view."""
        self._check()
        ptr = lib.LLKA_findAtomView(
            byref(self._c),
            *_find_atom_criteria(
                label_atom_id, label_comp_id, label_asym_id, label_seq_id, label_alt_id, pdbx_PDB_ins_code, pdbx_PDB_model_num
            ),
        )
        if not ptr:
            return None
        address = ctypes.cast(ptr, ctypes.c_void_p).value
        atoms = self._c.atoms
        return next(idx for idx in range(self._c.nAtoms) if ctypes.cast(atoms[idx], ctypes.c_void_p).value == address)

    def to_structure(self) -> Structure:
        """Returns a new structure with deep copies of the viewed atoms."""
        self._check()
        return Structure._take(lib.LLKA_makeStructureFromPtrs(self._c.atoms, self._c.nAtoms))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Structure | StructureView):
            return len(self) == len(other) and all(a == b for a, b in zip(self, other))
        return NotImplemented

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f'<StructureView of {len(self)} atoms>'


StructureLike = Structure | StructureView | Iterable[Atom]


def as_structure(stru: StructureLike) -> Structure:
    """Returns ``stru`` if it is a :class:`Structure`, otherwise creates a structure from the atoms."""
    if isinstance(stru, Structure):
        return stru
    if isinstance(stru, StructureView):
        return stru.to_structure()
    return Structure(stru)


def require_structure(stru: object, what: str) -> Structure:
    """Raises TypeError unless ``stru`` is a :class:`Structure`. Used where a view of ``stru`` is created."""
    if not isinstance(stru, Structure):
        raise TypeError(f'{what} requires a Structure, got {type(stru).__name__}')
    return stru


def checked_view(view: StructureView) -> LLKA_StructureView:
    view._check()
    return view._c


def base_kind(comp_id: str) -> BaseKind:
    """Returns whether a residue has a purine, pyrimidine or a non-standard base.

    :raises InvalidArgumentError: The residue name is not recognized.
    """
    kind = _capi.c_enum()
    code = lib.LLKA_baseKind(_capi.encode(comp_id), byref(kind))
    if code != 0:
        raise InvalidArgumentError(f'Unknown residue {comp_id!r}')
    return BaseKind(kind.value)


def compare_atoms(a: Atom, b: Atom, *, ignore_id: bool = False) -> bool:
    """Compares two atoms the way libLLKA does (``LLKA_compareAtoms``).

    Unlike ``a == b``, ``auth_`` identifiers are not compared and the ``id`` may be ignored.
    """
    ca = a._to_c()
    cb = b._to_c()
    return bool(lib.LLKA_compareAtoms(byref(ca), byref(cb), bool(ignore_id)))


