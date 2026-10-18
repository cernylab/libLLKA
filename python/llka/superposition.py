"""Centroids, RMSD and Kabsch superposition of point sets and structures (``llka_superposition.h``).

Functions that move coordinates (:func:`superpose`, :func:`apply_transformation`)
modify their first argument in place. Point sets passed as plain sequences are
therefore only accepted where they are not modified.
"""

from __future__ import annotations

from collections.abc import Iterable
from ctypes import byref, c_double

from . import _capi
from ._capi import LLKA_Matrix, lib
from .core import Matrix, Point, PointLike, Points, as_points
from .structure import Structure, StructureView, checked_view

__all__ = [
    'apply_transformation',
    'centroid',
    'rmsd',
    'superpose',
    'superposition_matrix',
]


def _check_transformation(matrix: Matrix) -> None:
    if not isinstance(matrix, Matrix):
        raise TypeError(f'Expected Matrix, got {type(matrix).__name__}')


def apply_transformation(what: Points | Structure, matrix: Matrix) -> None:
    """Applies a 4x4 transformation matrix to points or atoms of a structure in place.

    :raises InvalidArgumentError: The matrix is not a 4x4 matrix.
    """
    _check_transformation(matrix)
    if isinstance(what, Points):
        code = lib.LLKA_applyTransformationPoints(byref(what._c), byref(matrix._c))
    elif isinstance(what, Structure):
        code = lib.LLKA_applyTransformationStructure(byref(what._c), byref(matrix._c))
    else:
        raise TypeError(f'Expected Points or Structure, got {type(what).__name__}')
    _capi.check(code, 'Cannot apply transformation')


def centroid(what: Points | Structure | Iterable[PointLike]) -> Point:
    """Returns the centroid of a set of points or of the atoms of a structure."""
    if isinstance(what, Structure):
        return Point._from_c(lib.LLKA_centroidStructure(byref(what._c)))
    pts = as_points(what)
    return Point._from_c(lib.LLKA_centroidPoints(byref(pts._c)))


def rmsd(a: Points | Structure | Iterable[PointLike], b: Points | Structure | Iterable[PointLike]) -> float:
    """Returns the root mean square distance between two equally ordered sets of points or structures.

    :raises MismatchingSizesError: The sets do not have the same number of points.
    """
    result = c_double()
    if isinstance(a, Structure) and isinstance(b, Structure):
        code = lib.LLKA_rmsdStructures(byref(a._c), byref(b._c), byref(result))
    elif isinstance(a, Structure) or isinstance(b, Structure):
        raise TypeError('Either both or none of the arguments must be Structures')
    else:
        pa, pb = as_points(a), as_points(b)
        code = lib.LLKA_rmsdPoints(byref(pa._c), byref(pb._c), byref(result))
    _capi.check(code, 'Cannot calculate RMSD')
    return result.value


def superpose(what: Points | Structure, onto: Points | Structure | StructureView | Iterable[PointLike]) -> float:
    """Superposes ``what`` onto ``onto`` with the Kabsch algorithm. ``what`` is modified in place.

    Both arguments must have the same number of points (atoms) in the same order.

    :returns: RMSD of the superposed points.
    :raises MismatchingSizesError: The arguments do not have the same number of points.
    """
    result = c_double()
    if isinstance(what, Points):
        if isinstance(onto, Structure | StructureView):
            raise TypeError('Points can only be superposed onto points')
        onto_pts = as_points(onto)
        code = lib.LLKA_superposePoints(byref(what._c), byref(onto_pts._c), byref(result))
    elif isinstance(what, Structure):
        if isinstance(onto, Structure):
            code = lib.LLKA_superposeStructures(byref(what._c), byref(onto._c), byref(result))
        elif isinstance(onto, StructureView):
            code = lib.LLKA_superposeStructuresView(byref(what._c), byref(checked_view(onto)), byref(result))
        else:
            raise TypeError('A Structure can only be superposed onto a Structure or a StructureView')
    else:
        raise TypeError(f'Expected Points or Structure, got {type(what).__name__}')
    _capi.check(code, 'Cannot superpose')
    return result.value


def superposition_matrix(
    what: Points | Structure | StructureView | Iterable[PointLike],
    onto: Points | Structure | StructureView | Iterable[PointLike],
) -> Matrix:
    """Calculates the 4x4 transformation matrix that superposes ``what`` onto ``onto``.

    Both arguments must be of the same kind (points, structures or views) and have
    the same, non-zero number of points in the same order.

    :raises MismatchingSizesError: The arguments do not have the same number of points.
    :raises InvalidArgumentError: The arguments are empty.
    """
    matrix = LLKA_Matrix()
    if isinstance(what, Structure) and isinstance(onto, Structure):
        code = lib.LLKA_superpositionMatrixStructures(byref(what._c), byref(onto._c), byref(matrix))
    elif isinstance(what, StructureView) and isinstance(onto, StructureView):
        code = lib.LLKA_superpositionMatrixStructureViews(
            byref(checked_view(what)), byref(checked_view(onto)), byref(matrix)
        )
    elif isinstance(what, Structure | StructureView) or isinstance(onto, Structure | StructureView):
        raise TypeError('Both arguments must be Structures, StructureViews or point sets')
    else:
        pw, po = as_points(what), as_points(onto)
        code = lib.LLKA_superpositionMatrixPoints(byref(pw._c), byref(po._c), byref(matrix))
    _capi.check(code, 'Cannot calculate superposition matrix')
    return Matrix._from_c(matrix)

