"""Basic geometric types: points, point sets and matrices (``llka_main.h``)."""

from __future__ import annotations

import ctypes
import weakref
from collections.abc import Iterable, Iterator, MutableSequence, Sequence
from typing import NamedTuple, SupportsFloat, overload

from . import _capi
from ._capi import LLKA_Matrix, LLKA_Point, LLKA_Points, lib
from .enums import RetCode

__all__ = ['Matrix', 'Point', 'Points', 'error_to_string']


def error_to_string(retcode: RetCode | int) -> str:
    """Returns the name of a libLLKA return code."""
    return _capi.error_to_string(int(retcode))


class Point(NamedTuple):
    """Position in 3D space in cartesian coordinates."""

    x: float
    y: float
    z: float

    @classmethod
    def _from_c(cls, c: LLKA_Point) -> Point:
        return cls(c.x, c.y, c.z)

    def _to_c(self) -> LLKA_Point:
        return LLKA_Point(self.x, self.y, self.z)


PointLike = Point | tuple[SupportsFloat, SupportsFloat, SupportsFloat] | Sequence[SupportsFloat]


def _to_point(p: PointLike) -> Point:
    if isinstance(p, Point):
        return p
    x, y, z = p
    return Point(float(x), float(y), float(z))


class Points(MutableSequence[Point]):
    """Set of points in 3D space backed by memory managed by libLLKA (``LLKA_Points``).

    ``Points`` behaves like a fixed-size list of :class:`Point` objects.
    Superposition functions modify ``Points`` in place.
    """

    __slots__ = ('_c', '_finalizer', '__weakref__')

    def __init__(self, points: Iterable[PointLike] | int = ()) -> None:
        """Creates a set of points.

        :param points: Either an iterable of points or a number of points to create.
                       In the latter case, all coordinates are initialized to zero.
        """
        if isinstance(points, int):
            if points < 0:
                raise ValueError('Number of points must not be negative')
            c = lib.LLKA_zeroPoints(points)
        else:
            pts = [_to_point(p) for p in points]
            arr = (LLKA_Point * len(pts))(*[(p.x, p.y, p.z) for p in pts])
            c = lib.LLKA_makePoints(arr, len(pts))

        self._c = c
        self._finalizer = weakref.finalize(self, _destroy_points, c)

    def __len__(self) -> int:
        return self._c.nPoints

    def _index(self, idx: int) -> int:
        n = self._c.nPoints
        if idx < 0:
            idx += n
        if not 0 <= idx < n:
            raise IndexError('Points index out of range')
        return idx

    @overload
    def __getitem__(self, idx: int) -> Point: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Point]: ...

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return [self[i] for i in range(*idx.indices(len(self)))]
        return Point._from_c(self._c.points[self._index(idx)])

    def __setitem__(self, idx, value) -> None:
        if isinstance(idx, slice):
            indices = range(*idx.indices(len(self)))
            values = list(value)
            if len(values) != len(indices):
                raise ValueError('Points have fixed size, slice assignment must not change the number of points')
            for i, v in zip(indices, values):
                self[i] = v
            return

        p = _to_point(value)
        self._c.points[self._index(idx)] = LLKA_Point(p.x, p.y, p.z)

    def __delitem__(self, idx) -> None:
        raise TypeError('Points have fixed size')

    def insert(self, index: int, value: PointLike) -> None:
        raise TypeError('Points have fixed size')

    def __iter__(self) -> Iterator[Point]:
        c = self._c
        for idx in range(c.nPoints):
            yield Point._from_c(c.points[idx])

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Points | list | tuple):
            return list(self) == list(other)
        return NotImplemented

    def __repr__(self) -> str:
        return f'Points({list(self)!r})'

    def copy(self) -> Points:
        return Points(self)


def _destroy_points(c: LLKA_Points) -> None:
    if c.nPoints > 0:
        lib.LLKA_destroyPoints(ctypes.byref(c))


def as_points(points: Points | Iterable[PointLike]) -> Points:
    return points if isinstance(points, Points) else Points(points)


class Matrix:
    """Column-major matrix of doubles backed by memory managed by libLLKA (``LLKA_Matrix``).

    Elements are accessed as ``matrix[row, column]``.
    """

    __slots__ = ('_c', '_finalizer', '__weakref__')

    def __init__(self, rows: int | Iterable[Iterable[SupportsFloat]], columns: int | None = None) -> None:
        """Creates a matrix.

        :param rows: Either the number of rows of a zero-initialized matrix or a list of rows
                     with the values of the matrix.
        :param columns: Number of columns. Required if ``rows`` is a number.
        """
        if isinstance(rows, int):
            if columns is None:
                raise TypeError('Number of columns is required when the number of rows is given')
            values = None
            n_rows, n_cols = rows, columns
        else:
            if columns is not None:
                raise TypeError('Number of columns must not be given together with matrix values')
            values = [[float(v) for v in row] for row in rows]
            n_rows = len(values)
            n_cols = len(values[0]) if values else 0
            if any(len(row) != n_cols for row in values):
                raise ValueError('All rows of a matrix must have the same length')

        if n_rows < 1 or n_cols < 1:
            raise ValueError('Matrix must have at least one row and one column')

        c = LLKA_Matrix()
        lib.LLKA_initMatrix(n_rows, n_cols, ctypes.byref(c))
        self._adopt(c)

        for r in range(n_rows):
            for col in range(n_cols):
                self[r, col] = values[r][col] if values is not None else 0.0

    def _adopt(self, c: LLKA_Matrix) -> None:
        self._c = c
        self._finalizer = weakref.finalize(self, lib.LLKA_destroyMatrix, ctypes.byref(c))

    @classmethod
    def _from_c(cls, c: LLKA_Matrix) -> Matrix:
        """Takes ownership of a matrix created by libLLKA."""
        m = cls.__new__(cls)
        m._adopt(c)
        return m

    @property
    def rows(self) -> int:
        return self._c.nRows

    @property
    def columns(self) -> int:
        return self._c.nCols

    @property
    def shape(self) -> tuple[int, int]:
        return self._c.nRows, self._c.nCols

    def _offset(self, key: tuple[int, int]) -> int:
        row, col = key
        n_rows, n_cols = self._c.nRows, self._c.nCols
        if not (0 <= row < n_rows and 0 <= col < n_cols):
            raise IndexError(f'Matrix index {key} out of range for a {n_rows}x{n_cols} matrix')
        # Same as LLKA_matrixGet() and LLKA_matrixSet()
        return n_rows * col + row

    def __getitem__(self, key: tuple[int, int]) -> float:
        return self._c.data[self._offset(key)]

    def __setitem__(self, key: tuple[int, int], value: SupportsFloat) -> None:
        self._c.data[self._offset(key)] = float(value)

    def tolist(self) -> list[list[float]]:
        """Returns the matrix as a list of rows."""
        return [[self[r, c] for c in range(self.columns)] for r in range(self.rows)]

    def copy(self) -> Matrix:
        return Matrix._from_c(lib.LLKA_duplicateMatrix(ctypes.byref(self._c)))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Matrix):
            return self.tolist() == other.tolist()
        return NotImplemented

    def __repr__(self) -> str:
        return f'Matrix({self.tolist()!r})'
