"""Angles, dihedral angles and distances between atoms (``llka_measurements.h``)."""

from __future__ import annotations

from ctypes import byref

from ._capi import lib
from .structure import Atom

__all__ = ['measure_angle', 'measure_dihedral', 'measure_distance']


def measure_angle(a: Atom, b: Atom, c: Atom, *, single_precision: bool = False) -> float:
    """Returns the angle between atoms ``a``, ``b`` and ``c`` in radians, with ``b`` at the vertex.

    :param single_precision: Calculate in single precision (``LLKA_measureAnglef``).
    """
    ca, cb, cc = a._to_c(), b._to_c(), c._to_c()
    fn = lib.LLKA_measureAnglef if single_precision else lib.LLKA_measureAngle
    return fn(byref(ca), byref(cb), byref(cc))


def measure_dihedral(a: Atom, b: Atom, c: Atom, d: Atom, *, single_precision: bool = False) -> float:
    """Returns the dihedral angle defined by atoms ``a``, ``b``, ``c`` and ``d`` in radians in the range [-PI; PI].

    :param single_precision: Calculate in single precision (``LLKA_measureDihedralf``).
    """
    ca, cb, cc, cd = a._to_c(), b._to_c(), c._to_c(), d._to_c()
    fn = lib.LLKA_measureDihedralf if single_precision else lib.LLKA_measureDihedral
    return fn(byref(ca), byref(cb), byref(cc), byref(cd))


def measure_distance(a: Atom, b: Atom, *, single_precision: bool = False) -> float:
    """Returns the distance between atoms ``a`` and ``b`` in Ångströms.

    :param single_precision: Calculate in single precision (``LLKA_measureDistancef``).
    """
    ca, cb = a._to_c(), b._to_c()
    fn = lib.LLKA_measureDistancef if single_precision else lib.LLKA_measureDistance
    return fn(byref(ca), byref(cb))
