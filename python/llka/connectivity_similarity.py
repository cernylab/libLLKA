"""Connectivity and similarity of NtC steps (``llka_connectivity_similarity.h``)."""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Sequence
from ctypes import byref

from . import _capi
from ._capi import (
    LLKA_Connectivities,
    LLKA_Connectivity,
    LLKA_Similarities,
    LLKA_Similarity,
    LLKA_Structure,
    LLKA_Structures,
    c_enum,
    lib,
)
from .enums import NtC
from .ntc import _valid_ntc
from .structure import StructureLike, as_structure

__all__ = [
    'Connectivity',
    'Similarity',
    'measure_step_connectivity_ntcs',
    'measure_step_connectivity_ntcs_multiple_first',
    'measure_step_connectivity_ntcs_multiple_second',
    'measure_step_connectivity_structures',
    'measure_step_connectivity_structures_multiple',
    'measure_step_similarity_ntc',
    'measure_step_similarity_ntc_multiple',
    'measure_step_similarity_structure',
]


@dataclasses.dataclass(frozen=True, slots=True)
class Connectivity:
    """Connectivity of two consecutive steps.

    Distances between the C5' (O3') atoms of the second residue of the first
    step and the C5' (O3') atoms of the first residue of the second step.
    """

    C5_prime_distance: float
    O3_prime_distance: float

    @classmethod
    def _from_c(cls, c: LLKA_Connectivity) -> Connectivity:
        return cls(c.C5PrimeDistance, c.O3PrimeDistance)


@dataclasses.dataclass(frozen=True, slots=True)
class Similarity:
    """Similarity of two steps: RMSD after superposition and euclidean distance of the step metrics."""

    rmsd: float
    euclidean_distance: float

    @classmethod
    def _from_c(cls, c: LLKA_Similarity) -> Similarity:
        return cls(c.rmsd, c.euclideanDistance)


def _ntc_list(ntcs: Iterable[NtC]):
    """Returns a C array of NtCs terminated by NtC.INVALID and the number of NtCs."""
    valid = [_valid_ntc(ntc) for ntc in ntcs]
    return (c_enum * (len(valid) + 1))(*valid, NtC.INVALID), len(valid)


def _connectivities(n: int):
    arr = (LLKA_Connectivity * max(n, 1))()
    return arr, LLKA_Connectivities(arr, n)


def measure_step_connectivity_ntcs(
    position_first: StructureLike, ntc_first: NtC, position_second: StructureLike, ntc_second: NtC
) -> Connectivity:
    """Measures connectivity between two reference NtC steps superposed onto two real steps.

    :param position_first: Step that determines the position of the first reference step.
    :param ntc_first: NtC of the first reference step.
    :param position_second: Step that determines the position of the second reference step.
    :param ntc_second: NtC of the second reference step.
    """
    ntc_first, ntc_second = _valid_ntc(ntc_first), _valid_ntc(ntc_second)
    pf, ps = as_structure(position_first), as_structure(position_second)
    result = LLKA_Connectivity()
    _capi.check(
        lib.LLKA_measureStepConnectivityNtCs(byref(pf._c), ntc_first, byref(ps._c), ntc_second, byref(result)),
        'Cannot measure step connectivity',
    )
    return Connectivity._from_c(result)


def measure_step_connectivity_ntcs_multiple_first(
    position_first: StructureLike, ntcs_first: Iterable[NtC], position_second: StructureLike, ntc_second: NtC
) -> list[Connectivity]:
    """Like :func:`measure_step_connectivity_ntcs` for multiple NtCs of the first step.

    :returns: Connectivities in the order of ``ntcs_first``.
    """
    ntc_second = _valid_ntc(ntc_second)
    c_ntcs, n = _ntc_list(ntcs_first)
    if n == 0:
        return []
    pf, ps = as_structure(position_first), as_structure(position_second)
    arr, results = _connectivities(n)
    _capi.check(
        lib.LLKA_measureStepConnectivityNtCsMultipleFirst(byref(pf._c), c_ntcs, byref(ps._c), ntc_second, byref(results)),
        'Cannot measure step connectivity',
    )
    return [Connectivity._from_c(arr[idx]) for idx in range(n)]


def measure_step_connectivity_ntcs_multiple_second(
    position_first: StructureLike, ntc_first: NtC, position_second: StructureLike, ntcs_second: Iterable[NtC]
) -> list[Connectivity]:
    """Like :func:`measure_step_connectivity_ntcs` for multiple NtCs of the second step.

    :returns: Connectivities in the order of ``ntcs_second``.
    """
    ntc_first = _valid_ntc(ntc_first)
    c_ntcs, n = _ntc_list(ntcs_second)
    if n == 0:
        return []
    pf, ps = as_structure(position_first), as_structure(position_second)
    arr, results = _connectivities(n)
    _capi.check(
        lib.LLKA_measureStepConnectivityNtCsMultipleSecond(byref(pf._c), ntc_first, byref(ps._c), c_ntcs, byref(results)),
        'Cannot measure step connectivity',
    )
    return [Connectivity._from_c(arr[idx]) for idx in range(n)]


def measure_step_connectivity_structures(
    position_first: StructureLike, dinu_first: StructureLike, position_second: StructureLike, dinu_second: StructureLike
) -> Connectivity:
    """Measures connectivity between two dinucleotides placed as two consecutive steps.

    :param position_first: Step that determines the position of the first dinucleotide.
    :param dinu_first: Dinucleotide used as the first step, usually a reference NtC step.
    :param position_second: Step that determines the position of the second dinucleotide.
    :param dinu_second: Dinucleotide used as the second step, usually a reference NtC step.
    """
    pf, df = as_structure(position_first), as_structure(dinu_first)
    ps, ds = as_structure(position_second), as_structure(dinu_second)
    result = LLKA_Connectivity()
    _capi.check(
        lib.LLKA_measureStepConnectivityStructures(byref(pf._c), byref(df._c), byref(ps._c), byref(ds._c), byref(result)),
        'Cannot measure step connectivity',
    )
    return Connectivity._from_c(result)


def measure_step_connectivity_structures_multiple(
    position_first: StructureLike,
    dinu_first: StructureLike,
    position_second: StructureLike,
    dinus_second: Sequence[StructureLike],
) -> list[Connectivity]:
    """Like :func:`measure_step_connectivity_structures` for multiple second dinucleotides.

    :returns: Connectivities in the order of ``dinus_second``.
    """
    seconds = [as_structure(s) for s in dinus_second]
    if not seconds:
        return []
    pf, df, ps = as_structure(position_first), as_structure(dinu_first), as_structure(position_second)

    # The structures stay owned by their Python objects; the array only holds shallow copies of the structs
    c_seconds = (LLKA_Structure * len(seconds))(*[LLKA_Structure(s._c.atoms, s._c.nAtoms) for s in seconds])
    c_strus = LLKA_Structures(c_seconds, len(seconds))
    arr, results = _connectivities(len(seconds))
    _capi.check(
        lib.LLKA_measureStepConnectivityStructuresMultiple(
            byref(pf._c), byref(df._c), byref(ps._c), byref(c_strus), byref(results)
        ),
        'Cannot measure step connectivity',
    )
    return [Connectivity._from_c(arr[idx]) for idx in range(len(seconds))]


def measure_step_similarity_ntc(stru: StructureLike, ntc: NtC) -> Similarity:
    """Measures similarity between a step and a reference NtC."""
    ntc = _valid_ntc(ntc)
    stru = as_structure(stru)
    result = LLKA_Similarity()
    _capi.check(
        lib.LLKA_measureStepSimilarityNtC(byref(stru._c), ntc, byref(result)),
        'Cannot measure step similarity',
    )
    return Similarity._from_c(result)


def measure_step_similarity_ntc_multiple(stru: StructureLike, ntcs: Iterable[NtC]) -> list[Similarity]:
    """Measures similarity between a step and multiple reference NtCs.

    :returns: Similarities in the order of ``ntcs``.
    """
    c_ntcs, n = _ntc_list(ntcs)
    if n == 0:
        return []
    stru = as_structure(stru)
    arr = (LLKA_Similarity * n)()
    results = LLKA_Similarities(arr, n)
    _capi.check(
        lib.LLKA_measureStepSimilarityNtCMultiple(byref(stru._c), c_ntcs, byref(results)),
        'Cannot measure step similarity',
    )
    return [Similarity._from_c(arr[idx]) for idx in range(n)]


def measure_step_similarity_structure(stru: StructureLike, ref: StructureLike) -> Similarity:
    """Measures similarity between two steps. The metrics of ``ref`` are calculated from its structure."""
    stru, ref = as_structure(stru), as_structure(ref)
    result = LLKA_Similarity()
    _capi.check(
        lib.LLKA_measureStepSimilarityStructure(byref(stru._c), byref(ref._c), byref(result)),
        'Cannot measure step similarity',
    )
    return Similarity._from_c(result)

