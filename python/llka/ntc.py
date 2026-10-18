"""NtC steps: step detection, metrics, backbone extraction and NtC/CANA names (``llka_ntc.h``)."""

from __future__ import annotations

import dataclasses
from ctypes import byref
from typing import NamedTuple

from . import _capi
from ._capi import (
    LLKA_AtomNameQuad,
    LLKA_BackboneAtom,
    LLKA_StepInfo,
    LLKA_StepMetrics,
    LLKA_Structure,
    LLKA_StructureView,
    STEP_METRICS_FIELDS,
    lib,
)
from .enums import CANA, BaseKind, CrossResidueMetric, DinucleotideTorsion, NtC, ResidueInStep
from .errors import InvalidArgumentError, LLKAError
from .structure import Structure, StructureLike, StructureView, as_structure, require_structure

__all__ = [
    'INTERNAL_CANA_VERSION',
    'INTERNAL_NTC_VERSION',
    'AtomNameQuad',
    'BackboneAtom',
    'StepInfo',
    'StepMetrics',
    'backbone_atom_index',
    'calculate_step_metrics',
    'calculate_step_metrics_difference_against_reference',
    'cana_to_name',
    'cross_residue_metric',
    'cross_residue_metric_atoms_from_bases',
    'cross_residue_metric_atoms_from_structure',
    'cross_residue_metric_name',
    'dinucleotide_torsion',
    'dinucleotide_torsion_atoms_from_bases',
    'dinucleotide_torsion_atoms_from_structure',
    'dinucleotide_torsion_name',
    'extract_backbone',
    'extract_backbone_view',
    'extract_extended_backbone',
    'extract_extended_backbone_view',
    'extract_metrics_structure',
    'is_step',
    'name_to_cana',
    'name_to_ntc',
    'ntc_structure',
    'ntc_to_name',
    'structure_is_step',
]

INTERNAL_NTC_VERSION = '3.5'
INTERNAL_CANA_VERSION = '2.3'


class AtomNameQuad(NamedTuple):
    """Names of atoms that define a torsion or a distance. ``c`` and ``d`` are ``None`` for distances."""

    a: str
    b: str
    c: str | None
    d: str | None

    @classmethod
    def _from_c(cls, c: LLKA_AtomNameQuad) -> AtomNameQuad:
        def name(v: bytes) -> str | None:
            return v.decode('utf-8') if v else None

        return cls(c.a.decode('utf-8'), c.b.decode('utf-8'), name(c.c), name(c.d))


@dataclasses.dataclass(frozen=True, slots=True)
class BackboneAtom:
    """Identifies a backbone atom in a step."""

    residue: ResidueInStep
    name: str


@dataclasses.dataclass(frozen=True, slots=True)
class StepInfo:
    """Common properties of a NtC step."""

    first_seq_id: int
    first_seq_id_auth: int
    second_seq_id: int
    second_seq_id_auth: int
    first_base_kind: BaseKind
    second_base_kind: BaseKind

    @classmethod
    def _from_c(cls, c: LLKA_StepInfo) -> StepInfo:
        return cls(
            c.firstSeqId,
            c.firstSeqIdAuth,
            c.secondSeqId,
            c.secondSeqIdAuth,
            BaseKind(c.firstBaseKind),
            BaseKind(c.secondBaseKind),
        )


@dataclasses.dataclass(frozen=True, slots=True)
class StepMetrics:
    """Torsions and distances used to classify a step.

    Torsions are in radians in the range [-PI; PI], distances are in Ångströms.
    """

    delta_1: float
    epsilon_1: float
    zeta_1: float
    alpha_2: float
    beta_2: float
    gamma_2: float
    delta_2: float
    chi_1: float
    chi_2: float
    CC: float
    NN: float
    mu: float

    @classmethod
    def _from_c(cls, c: LLKA_StepMetrics) -> StepMetrics:
        return cls(*(getattr(c, name) for name in STEP_METRICS_FIELDS))

    def _to_c(self) -> LLKA_StepMetrics:
        return LLKA_StepMetrics(*(getattr(self, name) for name in STEP_METRICS_FIELDS))

    def get_torsion(self, torsion: DinucleotideTorsion) -> float:
        """Returns the value of a dinucleotide torsion."""
        return getattr(self, STEP_METRICS_FIELDS[DinucleotideTorsion(torsion)])


def _valid_ntc(ntc: NtC | int) -> NtC:
    try:
        ntc = NtC(ntc)
    except ValueError:
        raise InvalidArgumentError(f'{ntc!r} is not a valid NtC') from None
    if ntc == NtC.INVALID:
        raise InvalidArgumentError('NtC.INVALID is not a valid NtC')
    return ntc


def _view(fn, stru: Structure, context: str, func_name: str) -> StructureView:
    require_structure(stru, func_name)
    c = LLKA_StructureView()
    _capi.check(fn(byref(stru._c), byref(c)), context)
    return StructureView._take(c, stru)


def _extract(fn, stru: StructureLike, context: str) -> Structure:
    stru = as_structure(stru)
    c = LLKA_Structure()
    _capi.check(fn(byref(stru._c), byref(c)), context)
    return Structure._take(c)


def backbone_atom_index(atom: BackboneAtom, backbone: StructureLike) -> int:
    """Returns the index of a backbone atom in a structure that contains backbone atoms of a single step.

    :raises InvalidArgumentError: The atom is not present in the backbone.
    """
    backbone = as_structure(backbone)
    c_atom = LLKA_BackboneAtom(ResidueInStep(atom.residue), _capi.encode(atom.name))
    idx = lib.LLKA_backboneAtomIndex(byref(c_atom), byref(backbone._c))
    if idx == _capi.LLKA_INVALID_BKBN_ATOM_INDEX:
        raise InvalidArgumentError(f'Backbone atom {atom} was not found')
    return idx


def calculate_step_metrics(stru: StructureLike) -> StepMetrics:
    """Calculates the NtC metrics of a step. The structure must be a single valid step."""
    stru = as_structure(stru)
    c = LLKA_StepMetrics()
    _capi.check(lib.LLKA_calculateStepMetrics(byref(stru._c), byref(c)), 'Cannot calculate step metrics')
    return StepMetrics._from_c(c)


def calculate_step_metrics_difference_against_reference(stru: StructureLike, ntc: NtC) -> StepMetrics:
    """Calculates differences between the NtC metrics of a step and the metrics of a reference NtC.

    Torsion differences are in radians in the range [-PI; PI].
    """
    ntc = _valid_ntc(ntc)
    stru = as_structure(stru)
    c = LLKA_StepMetrics()
    _capi.check(
        lib.LLKA_calculateStepMetricsDifferenceAgainstReference(byref(stru._c), ntc, byref(c)),
        'Cannot calculate step metrics difference',
    )
    return StepMetrics._from_c(c)


def cana_to_name(cana: CANA) -> str:
    """Returns the human-readable name of a CANA class, e.g. ``"A-B"`` for ``CANA.AB``."""
    return lib.LLKA_CANAToName(CANA(cana)).decode('ascii')


def name_to_cana(name: str) -> CANA:
    """Converts a CANA name such as ``"BB2"`` to :class:`CANA`.

    :raises ValueError: The name is not a valid CANA name.
    """
    cana = CANA(lib.LLKA_nameToCANA(_capi.encode(name)))
    if cana == CANA.INVALID:
        raise ValueError(f'{name!r} is not a valid CANA name')
    return cana


def cross_residue_metric(metric: CrossResidueMetric, stru: StructureLike) -> Structure:
    """Returns copies of the atoms that define a cross-residue metric, in a fixed order.

    The structure must contain all atoms needed to calculate the step metrics,
    see :func:`extract_metrics_structure`.
    """
    stru = as_structure(stru)
    c = LLKA_Structure()
    _capi.check(
        lib.LLKA_crossResidueMetric(CrossResidueMetric(metric), byref(stru._c), byref(c)),
        'Cannot get cross-residue metric atoms',
    )
    return Structure._take(c)


def cross_residue_metric_atoms_from_bases(first_base: str, second_base: str, metric: CrossResidueMetric) -> AtomNameQuad:
    """Returns the names of atoms that define a cross-residue metric for a pair of residues."""
    quad = LLKA_AtomNameQuad()
    _capi.check(
        lib.LLKA_crossResidueMetricAtomsFromBases(
            _capi.encode(first_base), _capi.encode(second_base), CrossResidueMetric(metric), byref(quad)
        ),
        f'Cannot get cross-residue metric atoms for residues {first_base!r} and {second_base!r}',
    )
    return AtomNameQuad._from_c(quad)


def cross_residue_metric_atoms_from_structure(stru: StructureLike, metric: CrossResidueMetric) -> AtomNameQuad:
    """Returns the names of atoms that define a cross-residue metric for a step."""
    stru = as_structure(stru)
    quad = LLKA_AtomNameQuad()
    _capi.check(
        lib.LLKA_crossResidueMetricAtomsFromStructure(byref(stru._c), CrossResidueMetric(metric), byref(quad)),
        'Cannot get cross-residue metric atoms',
    )
    return AtomNameQuad._from_c(quad)


def cross_residue_metric_name(metric: CrossResidueMetric, greek: bool = False) -> str:
    """Returns the human-readable name of a cross-residue metric."""
    return lib.LLKA_crossResidueMetricName(CrossResidueMetric(metric), bool(greek)).decode('utf-8')


def dinucleotide_torsion(torsion: DinucleotideTorsion, stru: StructureLike) -> Structure:
    """Returns copies of the atoms that define a dinucleotide torsion, in a fixed order.

    The structure must contain all atoms needed to calculate the step metrics,
    see :func:`extract_metrics_structure`.
    """
    stru = as_structure(stru)
    c = LLKA_Structure()
    _capi.check(
        lib.LLKA_dinucleotideTorsion(DinucleotideTorsion(torsion), byref(stru._c), byref(c)),
        'Cannot get dinucleotide torsion atoms',
    )
    return Structure._take(c)


def dinucleotide_torsion_atoms_from_bases(
    first_base: str, second_base: str, torsion: DinucleotideTorsion
) -> AtomNameQuad:
    """Returns the names of atoms that define a dinucleotide torsion for a pair of residues."""
    quad = LLKA_AtomNameQuad()
    _capi.check(
        lib.LLKA_dinucleotideTorsionAtomsFromBases(
            _capi.encode(first_base), _capi.encode(second_base), DinucleotideTorsion(torsion), byref(quad)
        ),
        f'Cannot get dinucleotide torsion atoms for residues {first_base!r} and {second_base!r}',
    )
    return AtomNameQuad._from_c(quad)


def dinucleotide_torsion_atoms_from_structure(stru: StructureLike, torsion: DinucleotideTorsion) -> AtomNameQuad:
    """Returns the names of atoms that define a dinucleotide torsion for a step."""
    stru = as_structure(stru)
    quad = LLKA_AtomNameQuad()
    _capi.check(
        lib.LLKA_dinucleotideTorsionAtomsFromStructure(byref(stru._c), DinucleotideTorsion(torsion), byref(quad)),
        'Cannot get dinucleotide torsion atoms',
    )
    return AtomNameQuad._from_c(quad)


def dinucleotide_torsion_name(torsion: DinucleotideTorsion, greek: bool = False) -> str:
    """Returns the human-readable name of a dinucleotide torsion."""
    return lib.LLKA_dinucleotideTorsionName(DinucleotideTorsion(torsion), bool(greek)).decode('utf-8')


def extract_backbone(stru: StructureLike) -> Structure:
    """Returns copies of the backbone atoms of a step, in a fixed order."""
    return _extract(lib.LLKA_extractBackbone, stru, 'Cannot extract backbone')


def extract_backbone_view(stru: Structure) -> StructureView:
    """Like :func:`extract_backbone` but returns a view of the atoms of ``stru``."""
    return _view(lib.LLKA_extractBackboneView, stru, 'Cannot extract backbone', 'extract_backbone_view()')


def extract_extended_backbone(stru: StructureLike) -> Structure:
    """Returns copies of the extended backbone atoms of a step, in a fixed order.

    The extended backbone also contains the C1' atoms and two base atoms of each residue.
    """
    return _extract(lib.LLKA_extractExtendedBackbone, stru, 'Cannot extract extended backbone')


def extract_extended_backbone_view(stru: Structure) -> StructureView:
    """Like :func:`extract_extended_backbone` but returns a view of the atoms of ``stru``."""
    return _view(lib.LLKA_extractExtendedBackboneView, stru, 'Cannot extract extended backbone', 'extract_extended_backbone_view()')


def extract_metrics_structure(stru: StructureLike) -> Structure:
    """Returns copies of the atoms of a step that are needed to calculate the step metrics."""
    return _extract(lib.LLKA_extractMetricsStructure, stru, 'Cannot extract metrics structure')


def name_to_ntc(name: str) -> NtC:
    """Converts a NtC name such as ``"AA00"`` to :class:`NtC`.

    :raises ValueError: The name is not a valid NtC name.
    """
    ntc = NtC(lib.LLKA_nameToNtC(_capi.encode(name)))
    if ntc == NtC.INVALID:
        raise ValueError(f'{name!r} is not a valid NtC name')
    return ntc


def ntc_structure(ntc: NtC) -> Structure:
    """Returns the structure of the reference dinucleotide of a NtC."""
    return Structure._take(lib.LLKA_NtCStructure(_valid_ntc(ntc)))


def ntc_to_name(ntc: NtC) -> str:
    """Returns the name of a NtC, e.g. ``"AA00"``."""
    return lib.LLKA_NtCToName(NtC(ntc)).decode('ascii')


def structure_is_step(stru: StructureLike) -> StepInfo:
    """Checks that a structure is a single NtC step and returns information about the step.

    :raises LLKAError: The structure is not a valid step. The exception type describes the reason.
    """
    stru = as_structure(stru)
    info = LLKA_StepInfo()
    _capi.check(lib.LLKA_structureIsStep(byref(stru._c), byref(info)), 'Structure is not a valid step')
    return StepInfo._from_c(info)


def is_step(stru: StructureLike) -> bool:
    """Returns whether a structure is a single valid NtC step."""
    try:
        structure_is_step(stru)
    except LLKAError:
        return False
    return True
