"""Loading of classification parameters (``llka_resource_loaders.h``).

The loaded parameters are returned as a :class:`ResourceList`, a read-only
sequence that keeps the data in memory managed by libLLKA. A ``ResourceList``
can be passed to :class:`~llka.classification.ClassificationContext` without
any conversion.
"""

from __future__ import annotations

import ctypes
import os
import weakref
from collections.abc import Callable, Iterator, Sequence
from ctypes import POINTER, byref
from typing import Generic, TypeVar, overload

from . import _capi
from ._capi import (
    LLKA_ClassificationCluster,
    LLKA_ClusterNuAngles,
    LLKA_Confal,
    LLKA_ConfalPercentile,
    LLKA_GoldenStep,
    LLKA_Resource,
    lib,
)
from .classification import ClassificationCluster, ClusterNuAngles, Confal, GoldenStep
from .enums import ResourceType

__all__ = ['ResourceList', 'load_resource_file', 'load_resource_text']

T = TypeVar('T')

_RESOURCE_TYPES: dict[ResourceType, tuple[type, Callable]] = {
    ResourceType.AVERAGE_NU_ANGLES: (LLKA_ClusterNuAngles, ClusterNuAngles._from_c),
    ResourceType.CLUSTERS: (LLKA_ClassificationCluster, ClassificationCluster._from_c),
    ResourceType.CONFALS: (LLKA_Confal, Confal._from_c),
    ResourceType.GOLDEN_STEPS: (LLKA_GoldenStep, GoldenStep._from_c),
    ResourceType.CONFAL_PERCENTILES: (LLKA_ConfalPercentile, lambda c: c.value),
}


def _destroy_resource(c: LLKA_Resource) -> None:
    lib.LLKA_destroyResource(byref(c))


class ResourceList(Sequence[T], Generic[T]):
    """Read-only sequence of loaded classification parameters.

    Items are :class:`~llka.classification.ClusterNuAngles`,
    :class:`~llka.classification.ClassificationCluster`,
    :class:`~llka.classification.Confal`, :class:`~llka.classification.GoldenStep`
    or ``float`` (confal percentiles), depending on :attr:`resource_type`.
    """

    __slots__ = ('_c', '_c_type', '_convert', '_finalizer', 'resource_type', '__weakref__')

    def __init__(self) -> None:
        raise TypeError('ResourceList cannot be created directly, use load_resource_file() or load_resource_text()')

    @classmethod
    def _take(cls, c: LLKA_Resource) -> ResourceList:
        self = cls.__new__(cls)
        self.resource_type = ResourceType(c.type)
        self._c = c
        self._c_type, self._convert = _RESOURCE_TYPES[self.resource_type]
        self._finalizer = weakref.finalize(self, _destroy_resource, c)
        return self

    def _c_pointer(self):
        return ctypes.cast(self._c.data, POINTER(self._c_type))

    def __len__(self) -> int:
        return self._c.count

    @overload
    def __getitem__(self, idx: int) -> T: ...
    @overload
    def __getitem__(self, idx: slice) -> list[T]: ...

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return [self[i] for i in range(*idx.indices(len(self)))]
        n = len(self)
        if idx < 0:
            idx += n
        if not 0 <= idx < n:
            raise IndexError('ResourceList index out of range')
        return self._convert(self._c_pointer()[idx])

    def __iter__(self) -> Iterator[T]:
        ptr = self._c_pointer()
        for idx in range(len(self)):
            yield self._convert(ptr[idx])

    def __repr__(self) -> str:
        return f'<ResourceList of {len(self)} {self.resource_type.name}>'


def load_resource_file(path: str | os.PathLike[str], resource_type: ResourceType) -> ResourceList:
    """Loads classification parameters from a file.

    :raises NoFileError: The file does not exist.
    :raises CannotReadFileError: The file cannot be read.
    :raises BadDataError: The file contains malformed data.
    """
    fn = _capi.require_function('LLKA_loadResourceFile')
    c = LLKA_Resource()
    c.type = ResourceType(resource_type)
    _capi.check(fn(_capi.path_arg(path), byref(c)), f'Cannot load {ResourceType(resource_type).name} from {os.fspath(path)}')
    return ResourceList._take(c)


def load_resource_text(text: str, resource_type: ResourceType) -> ResourceList:
    """Loads classification parameters from a string with the contents of a parameters file.

    :raises BadDataError: The text contains malformed data.
    """
    c = LLKA_Resource()
    c.type = ResourceType(resource_type)
    _capi.check(
        lib.LLKA_loadResourceText(_capi.encode(text), byref(c)),
        f'Cannot load {ResourceType(resource_type).name}',
    )
    return ResourceList._take(c)
