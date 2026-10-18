"""Nucleotides, ribose geometry and sugar puckers (``llka_nucleotide.h``)."""

from __future__ import annotations

import dataclasses
from ctypes import byref

from . import _capi
from ._capi import NU_ANGLES_FIELDS, LLKA_NuAngles, LLKA_RiboseMetrics, LLKA_Structure, LLKA_StructureView, lib
from .enums import SugarPucker, SugarPuckerNameBrevity
from .errors import InvalidArgumentError
from .structure import Structure, StructureLike, StructureView, as_structure, checked_view, require_structure

__all__ = [
    'NuAngles',
    'RiboseMetrics',
    'extract_nucleotide',
    'extract_nucleotide_view',
    'extract_ribose',
    'extract_ribose_view',
    'is_nucleotide_compound',
    'name_to_sugar_pucker',
    'ribose_metrics',
    'sugar_pucker',
    'sugar_pucker_to_name',
]


@dataclasses.dataclass(frozen=True, slots=True)
class NuAngles:
    """Torsion angles of the (deoxy)ribose ring in radians."""

    nu_0: float  # C4'-O4'-C1'-C2'
    nu_1: float  # O4'-C1'-C2'-C3'
    nu_2: float  # C1'-C2'-C3'-C4'
    nu_3: float  # C2'-C3'-C4'-O4'
    nu_4: float  # C3'-C4'-O4'-C1'

    @classmethod
    def _from_c(cls, c: LLKA_NuAngles) -> NuAngles:
        return cls(*(getattr(c, name) for name in NU_ANGLES_FIELDS))

    def _to_c(self) -> LLKA_NuAngles:
        return LLKA_NuAngles(*(getattr(self, name) for name in NU_ANGLES_FIELDS))


@dataclasses.dataclass(frozen=True, slots=True)
class RiboseMetrics:
    """Geometry of a (deoxy)ribose ring."""

    nus: NuAngles
    P: float          # Pseudorotation
    t_max: float      # tauMax
    pucker: SugarPucker

    @classmethod
    def _from_c(cls, c: LLKA_RiboseMetrics) -> RiboseMetrics:
        return cls(NuAngles._from_c(c.nus), c.P, c.tMax, SugarPucker(c.pucker))


def extract_nucleotide(stru: StructureLike, model_num: int, asym_id: str, seq_id: int) -> Structure:
    """Returns copies of the atoms of a single nucleotide. The result is empty if there is no such nucleotide.

    :param model_num: ``pdbx_PDB_model_num`` of the nucleotide.
    :param asym_id: ``label_asym_id`` of the nucleotide.
    :param seq_id: ``label_seq_id`` of the nucleotide.
    """
    stru = as_structure(stru)
    return Structure._take(lib.LLKA_extractNucleotide(byref(stru._c), model_num, _capi.encode(asym_id), seq_id))


def extract_nucleotide_view(stru: Structure, model_num: int, asym_id: str, seq_id: int) -> StructureView:
    """Like :func:`extract_nucleotide` but returns a view of the atoms of ``stru``."""
    require_structure(stru, 'extract_nucleotide_view()')
    c = lib.LLKA_extractNucleotideView(byref(stru._c), model_num, _capi.encode(asym_id), seq_id)
    return StructureView._take(c, stru)


def extract_ribose(stru: StructureLike) -> Structure:
    """Returns copies of the ribose ring atoms of a structure that contains exactly one ribose."""
    stru = as_structure(stru)
    c = LLKA_Structure()
    _capi.check(lib.LLKA_extractRibose(byref(stru._c), byref(c)), 'Cannot extract ribose')
    return Structure._take(c)


def extract_ribose_view(stru: Structure) -> StructureView:
    """Like :func:`extract_ribose` but returns a view of the atoms of ``stru``."""
    require_structure(stru, 'extract_ribose_view()')
    c = LLKA_StructureView()
    _capi.check(lib.LLKA_extractRiboseView(byref(stru._c), byref(c)), 'Cannot extract ribose')
    return StructureView._take(c, stru)


def is_nucleotide_compound(comp_id: str) -> bool:
    """Returns whether a compound name is recognized as a nucleotide."""
    return bool(lib.LLKA_isNucleotideCompound(_capi.encode(comp_id)))


def name_to_sugar_pucker(name: str) -> SugarPucker:
    """Converts a sugar pucker name such as ``"C1endo"`` or ``"C3' exo"`` to :class:`SugarPucker`.

    :raises ValueError: The name is not recognized.
    """
    pucker = SugarPucker(lib.LLKA_nameToSugarPucker(_capi.encode(name)))
    if pucker == SugarPucker.INVALID:
        raise ValueError(f'{name!r} is not a valid sugar pucker name')
    return pucker


def ribose_metrics(stru: StructureLike) -> RiboseMetrics:
    """Calculates the geometry of the ribose ring of a structure or view that contains exactly one ribose."""
    c = LLKA_RiboseMetrics()
    if isinstance(stru, StructureView):
        if len(stru) == 0:
            raise InvalidArgumentError('Cannot calculate ribose metrics of an empty view')
        _capi.check(lib.LLKA_riboseMetricsView(byref(checked_view(stru)), byref(c)), 'Cannot calculate ribose metrics')
    else:
        stru = as_structure(stru)
        if len(stru) == 0:
            raise InvalidArgumentError('Cannot calculate ribose metrics of an empty structure')
        _capi.check(lib.LLKA_riboseMetrics(byref(stru._c), byref(c)), 'Cannot calculate ribose metrics')
    return RiboseMetrics._from_c(c)


def sugar_pucker(stru: Structure | StructureView) -> SugarPucker:
    """Determines the sugar pucker of a structure that contains exactly one ribose."""
    pucker = _capi.c_enum()
    if isinstance(stru, StructureView):
        if len(stru) == 0:
            raise InvalidArgumentError('Cannot determine sugar pucker of an empty view')
        _capi.check(lib.LLKA_sugarPuckerView(byref(checked_view(stru)), byref(pucker)), 'Cannot determine sugar pucker')
    else:
        stru = as_structure(stru)
        if len(stru) == 0:
            raise InvalidArgumentError('Cannot determine sugar pucker of an empty structure')
        _capi.check(lib.LLKA_sugarPucker(byref(stru._c), byref(pucker)), 'Cannot determine sugar pucker')
    return SugarPucker(pucker.value)


def sugar_pucker_to_name(pucker: SugarPucker, brevity: SugarPuckerNameBrevity = SugarPuckerNameBrevity.TERSE) -> str:
    """Returns the name of a sugar pucker. ``SugarPucker.INVALID`` yields an empty string."""
    return lib.LLKA_sugarPuckerToName(SugarPucker(pucker), SugarPuckerNameBrevity(brevity)).decode('utf-8')
