"""NtC and CANA classification of dinucleotide steps (``llka_classification.h``).

A :class:`ClassificationContext` is created from the classification parameters,
usually loaded with :func:`llka.resources.load_resource_file`, and then used to
classify steps::

    ctx = ClassificationContext(
        clusters=load_resource_file('clusters.csv', ResourceType.CLUSTERS),
        golden_steps=load_resource_file('golden_steps.csv', ResourceType.GOLDEN_STEPS),
        confals=load_resource_file('confals.csv', ResourceType.CONFALS),
        cluster_nu_angles=load_resource_file('nu_angles.csv', ResourceType.AVERAGE_NU_ANGLES),
        confal_percentiles=load_resource_file('confal_percentiles.csv', ResourceType.CONFAL_PERCENTILES),
        limits=ClassificationLimits(...),
        max_close_enough_rmsd=0.5,
    )
    for step in structure.split_to_dinucleotide_steps():
        print(ctx.classify_step(step).assigned_ntc)
"""

from __future__ import annotations

import ctypes
import dataclasses
import math
import weakref
from collections.abc import Iterable, Sequence
from ctypes import POINTER, byref

from . import _capi
from ._capi import (
    NU_ANGLES_FIELDS,
    STEP_METRICS_FIELDS,
    LLKA_AttemptedClassifiedStep,
    LLKA_ClassificationCluster,
    LLKA_ClassificationContext,
    LLKA_ClassificationLimits,
    LLKA_ClassificationMetric,
    LLKA_ClassifiedStep,
    LLKA_ClassifiedSteps,
    LLKA_ClusterNuAngles,
    LLKA_Confal,
    LLKA_ConfalPercentile,
    LLKA_ConfalScore,
    LLKA_GoldenStep,
    LLKA_NuAnglesMetrics,
    LLKA_Structure,
    LLKA_Structures,
    lib,
)
from .enums import CANA, ClassificationViolation, NtC, RetCode, SugarPucker, ViolatingTorsion
from .errors import LLKAError, error_for_retcode
from .ntc import StepMetrics
from .nucleotide import NuAngles
from .structure import StructureLike, as_structure

__all__ = [
    'AttemptedClassifiedStep',
    'AverageConfal',
    'ClassificationCluster',
    'ClassificationContext',
    'ClassificationLimits',
    'ClassificationMetric',
    'ClassifiedStep',
    'ClusterNuAngles',
    'Confal',
    'ConfalScore',
    'GoldenStep',
    'NuAnglesMetrics',
    'classification_violation_to_name',
]


def _optional_ntc(value: int) -> NtC | None:
    ntc = NtC(value)
    return None if ntc == NtC.INVALID else ntc


def _optional_cana(value: int) -> CANA | None:
    cana = CANA(value)
    return None if cana == CANA.INVALID else cana


@dataclasses.dataclass(frozen=True, slots=True)
class ClassificationMetric:
    """Tolerances of a classification metric. Angles are in radians in the range [0; 2PI]."""

    deviation: float
    min_value: float
    mean_value: float
    max_value: float

    @classmethod
    def _from_c(cls, c: LLKA_ClassificationMetric) -> ClassificationMetric:
        return cls(c.deviation, c.minValue, c.meanValue, c.maxValue)

    def _to_c(self) -> LLKA_ClassificationMetric:
        return LLKA_ClassificationMetric(self.deviation, self.min_value, self.mean_value, self.max_value)


@dataclasses.dataclass(frozen=True, slots=True)
class NuAnglesMetrics:
    """Classification metrics of the (deoxy)ribose torsion angles."""

    nu_0: ClassificationMetric
    nu_1: ClassificationMetric
    nu_2: ClassificationMetric
    nu_3: ClassificationMetric
    nu_4: ClassificationMetric

    @classmethod
    def _from_c(cls, c: LLKA_NuAnglesMetrics) -> NuAnglesMetrics:
        return cls(*(ClassificationMetric._from_c(getattr(c, name)) for name in NU_ANGLES_FIELDS))

    def _to_c(self) -> LLKA_NuAnglesMetrics:
        return LLKA_NuAnglesMetrics(*(getattr(self, name)._to_c() for name in NU_ANGLES_FIELDS))


@dataclasses.dataclass(frozen=True, slots=True)
class ClusterNuAngles:
    """Nu angles metrics of both nucleotides of a classification cluster."""

    first_nucleotide: NuAnglesMetrics
    second_nucleotide: NuAnglesMetrics
    cluster_number: int

    @classmethod
    def _from_c(cls, c: LLKA_ClusterNuAngles) -> ClusterNuAngles:
        return cls(
            NuAnglesMetrics._from_c(c.firstNucleotide),
            NuAnglesMetrics._from_c(c.secondNucleotide),
            c.clusterNumber,
        )

    def _fill_c(self, c: LLKA_ClusterNuAngles) -> None:
        c.firstNucleotide = self.first_nucleotide._to_c()
        c.secondNucleotide = self.second_nucleotide._to_c()
        c.clusterNumber = self.cluster_number


@dataclasses.dataclass(frozen=True, slots=True)
class ClassificationCluster:
    """Classification cluster that corresponds to one NtC and CANA class.

    Torsion metrics are in radians, distance metrics in Ångströms.
    """

    delta_1: ClassificationMetric
    epsilon_1: ClassificationMetric
    zeta_1: ClassificationMetric
    alpha_2: ClassificationMetric
    beta_2: ClassificationMetric
    gamma_2: ClassificationMetric
    delta_2: ClassificationMetric
    chi_1: ClassificationMetric
    chi_2: ClassificationMetric
    CC: ClassificationMetric
    NN: ClassificationMetric
    mu: ClassificationMetric
    nus_first: NuAnglesMetrics
    nus_second: NuAnglesMetrics
    ribose_pseudorotation_1: float
    ribose_pseudorotation_2: float
    number: int
    ntc: NtC
    cana: CANA

    @classmethod
    def _from_c(cls, c: LLKA_ClassificationCluster) -> ClassificationCluster:
        return cls(
            *(ClassificationMetric._from_c(getattr(c, name)) for name in STEP_METRICS_FIELDS),
            nus_first=NuAnglesMetrics._from_c(c.nusFirst),
            nus_second=NuAnglesMetrics._from_c(c.nusSecond),
            ribose_pseudorotation_1=c.ribosePseudorotation_1,
            ribose_pseudorotation_2=c.ribosePseudorotation_2,
            number=c.number,
            ntc=NtC(c.NtC),
            cana=CANA(c.CANA),
        )

    def _fill_c(self, c: LLKA_ClassificationCluster) -> None:
        for name in STEP_METRICS_FIELDS:
            setattr(c, name, getattr(self, name)._to_c())
        c.nusFirst = self.nus_first._to_c()
        c.nusSecond = self.nus_second._to_c()
        c.ribosePseudorotation_1 = self.ribose_pseudorotation_1
        c.ribosePseudorotation_2 = self.ribose_pseudorotation_2
        c.number = self.number
        c.NtC = NtC(self.ntc)
        c.CANA = CANA(self.cana)


@dataclasses.dataclass(frozen=True, slots=True)
class Confal:
    """Confal definition of a classification cluster."""

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
    nus_first: NuAngles
    nus_second: NuAngles
    cluster_number: int

    @classmethod
    def _from_c(cls, c: LLKA_Confal) -> Confal:
        return cls(
            *(getattr(c, name) for name in STEP_METRICS_FIELDS),
            nus_first=NuAngles._from_c(c.nusFirst),
            nus_second=NuAngles._from_c(c.nusSecond),
            cluster_number=c.clusterNumber,
        )

    def _fill_c(self, c: LLKA_Confal) -> None:
        for name in STEP_METRICS_FIELDS:
            setattr(c, name, getattr(self, name))
        c.nusFirst = self.nus_first._to_c()
        c.nusSecond = self.nus_second._to_c()
        c.clusterNumber = self.cluster_number


@dataclasses.dataclass(frozen=True, slots=True)
class GoldenStep:
    """Golden step definition."""

    pucker_1: SugarPucker
    pucker_2: SugarPucker
    nu_angles_1: NuAngles
    nu_angles_2: NuAngles
    metrics: StepMetrics
    name: str
    cluster_number: int

    @classmethod
    def _from_c(cls, c: LLKA_GoldenStep) -> GoldenStep:
        return cls(
            SugarPucker(c.pucker_1),
            SugarPucker(c.pucker_2),
            NuAngles._from_c(c.nuAngles_1),
            NuAngles._from_c(c.nuAngles_2),
            StepMetrics._from_c(c.metrics),
            _capi.decode(c.name),
            c.clusterNumber,
        )

    def _fill_c(self, c: LLKA_GoldenStep) -> None:
        c.pucker_1 = SugarPucker(self.pucker_1)
        c.pucker_2 = SugarPucker(self.pucker_2)
        c.nuAngles_1 = self.nu_angles_1._to_c()
        c.nuAngles_2 = self.nu_angles_2._to_c()
        c.metrics = self.metrics._to_c()
        c.name = _capi.encode(self.name)
        c.clusterIdx = 0
        c.clusterNumber = self.cluster_number


@dataclasses.dataclass(frozen=True, slots=True)
class ClassificationLimits:
    """Tolerances and options of the classification. Angles are in radians."""

    minimum_nearest_neighbors: int
    number_of_used_nearest_neighbors: int
    average_neighbors_torsion_cutoff: float
    nearest_neighbor_torsions_cutoff: float
    total_distance_cutoff: float
    pseudorotation_cutoff: float
    minimum_cluster_votes: float

    def _to_c(self) -> LLKA_ClassificationLimits:
        return LLKA_ClassificationLimits(
            self.minimum_nearest_neighbors,
            self.number_of_used_nearest_neighbors,
            self.average_neighbors_torsion_cutoff,
            self.nearest_neighbor_torsions_cutoff,
            self.total_distance_cutoff,
            self.pseudorotation_cutoff,
            self.minimum_cluster_votes,
        )


@dataclasses.dataclass(frozen=True, slots=True)
class ConfalScore:
    """Confal scores of a classified step."""

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
    total: float

    @classmethod
    def _from_c(cls, c: LLKA_ConfalScore) -> ConfalScore:
        return cls(*(getattr(c, name) for name in STEP_METRICS_FIELDS), total=c.total)

    def _to_c(self) -> LLKA_ConfalScore:
        return LLKA_ConfalScore(*(getattr(self, name) for name in STEP_METRICS_FIELDS), self.total)


@dataclasses.dataclass(frozen=True, slots=True)
class AverageConfal:
    """Average confal score of multiple steps and its statistical percentile."""

    score: float
    percentile: float


@dataclasses.dataclass(frozen=True, slots=True)
class ClassifiedStep:
    """Result of an attempt to classify a step.

    ``assigned_ntc`` and ``assigned_cana`` are ``None`` unless the step was
    fully classified, i.e. unless ``violations`` is ``ClassificationViolation.OK``.
    """

    assigned_ntc: NtC | None
    assigned_cana: CANA | None
    closest_ntc: NtC | None
    closest_cana: CANA | None
    confal_score: ConfalScore
    euclidean_distance_ntc_ideal: float
    metrics: StepMetrics
    differences_from_ntc_averages: StepMetrics
    nu_angles_1: NuAngles
    nu_angles_2: NuAngles
    ribose_pseudorotation_1: float
    ribose_pseudorotation_2: float
    tau_1: float
    tau_2: float
    sugar_pucker_1: SugarPucker
    sugar_pucker_2: SugarPucker
    nu_angle_differences_1: NuAngles
    nu_angle_differences_2: NuAngles
    rmsd_to_closest_ntc: float
    closest_golden_step: str
    violations: ClassificationViolation
    violating_torsions_average: ViolatingTorsion
    violating_torsions_nearest: ViolatingTorsion

    @property
    def is_assigned(self) -> bool:
        """Whether NtC and CANA classes were assigned to the step."""
        return self.violations == ClassificationViolation.OK

    @classmethod
    def _from_c(cls, c: LLKA_ClassifiedStep) -> ClassifiedStep:
        return cls(
            assigned_ntc=_optional_ntc(c.assignedNtC),
            assigned_cana=_optional_cana(c.assignedCANA),
            closest_ntc=_optional_ntc(c.closestNtC),
            closest_cana=_optional_cana(c.closestCANA),
            confal_score=ConfalScore._from_c(c.confalScore),
            euclidean_distance_ntc_ideal=c.euclideanDistanceNtCIdeal,
            metrics=StepMetrics._from_c(c.metrics),
            differences_from_ntc_averages=StepMetrics._from_c(c.differencesFromNtCAverages),
            nu_angles_1=NuAngles._from_c(c.nuAngles_1),
            nu_angles_2=NuAngles._from_c(c.nuAngles_2),
            ribose_pseudorotation_1=c.ribosePseudorotation_1,
            ribose_pseudorotation_2=c.ribosePseudorotation_2,
            tau_1=c.tau_1,
            tau_2=c.tau_2,
            sugar_pucker_1=SugarPucker(c.sugarPucker_1),
            sugar_pucker_2=SugarPucker(c.sugarPucker_2),
            nu_angle_differences_1=NuAngles._from_c(c.nuAngleDifferences_1),
            nu_angle_differences_2=NuAngles._from_c(c.nuAngleDifferences_2),
            rmsd_to_closest_ntc=c.rmsdToClosestNtC,
            closest_golden_step=_capi.decode(c.closestGoldenStep),
            violations=ClassificationViolation(c.violations & 0xFFFFFFFF),
            violating_torsions_average=ViolatingTorsion(c.violatingTorsionsAverage & 0xFFFF),
            violating_torsions_nearest=ViolatingTorsion(c.violatingTorsionsNearest & 0xFFFF),
        )

    def _fill_c(self, c: LLKA_ClassifiedStep) -> None:
        def ntc(v: NtC | None) -> int:
            return NtC.INVALID if v is None else NtC(v)

        def cana(v: CANA | None) -> int:
            return CANA.INVALID if v is None else CANA(v)

        c.assignedNtC = ntc(self.assigned_ntc)
        c.assignedCANA = cana(self.assigned_cana)
        c.closestNtC = ntc(self.closest_ntc)
        c.closestCANA = cana(self.closest_cana)
        c.confalScore = self.confal_score._to_c()
        c.euclideanDistanceNtCIdeal = self.euclidean_distance_ntc_ideal
        c.metrics = self.metrics._to_c()
        c.differencesFromNtCAverages = self.differences_from_ntc_averages._to_c()
        c.nuAngles_1 = self.nu_angles_1._to_c()
        c.nuAngles_2 = self.nu_angles_2._to_c()
        c.ribosePseudorotation_1 = self.ribose_pseudorotation_1
        c.ribosePseudorotation_2 = self.ribose_pseudorotation_2
        c.tau_1 = self.tau_1
        c.tau_2 = self.tau_2
        c.sugarPucker_1 = SugarPucker(self.sugar_pucker_1)
        c.sugarPucker_2 = SugarPucker(self.sugar_pucker_2)
        c.nuAngleDifferences_1 = self.nu_angle_differences_1._to_c()
        c.nuAngleDifferences_2 = self.nu_angle_differences_2._to_c()
        c.rmsdToClosestNtC = self.rmsd_to_closest_ntc
        c.closestGoldenStep = _capi.encode(self.closest_golden_step)
        c.violations = ctypes.c_int32(int(self.violations)).value
        c.violatingTorsionsAverage = ctypes.c_int16(int(self.violating_torsions_average)).value
        c.violatingTorsionsNearest = ctypes.c_int16(int(self.violating_torsions_nearest)).value


@dataclasses.dataclass(frozen=True, slots=True)
class AttemptedClassifiedStep:
    """Result of an attempt to classify one of multiple steps.

    If the step could not be classified at all (e.g. because it was not a valid
    step), ``step`` is ``None`` and ``error`` holds the corresponding exception.
    """

    step: ClassifiedStep | None
    error: LLKAError | None

    @property
    def ok(self) -> bool:
        return self.error is None


def classification_violation_to_name(violation: ClassificationViolation | int) -> str:
    """Returns the name of a single classification violation flag."""
    return lib.LLKA_classificationViolationToName(ctypes.c_int32(int(violation)).value).decode('ascii')


def _to_c_array(items, c_type, context: str):
    """Returns a pointer to a C array of the items and the number of items.

    Resources loaded by :mod:`llka.resources` are passed to libLLKA without conversion.
    """
    from .resources import ResourceList

    if isinstance(items, ResourceList):
        if items._c_type is not c_type:
            raise TypeError(f'{context}: resource of a wrong type')
        return items._c_pointer(), len(items), items

    items = list(items)
    arr = (c_type * len(items))()
    for idx, item in enumerate(items):
        if c_type is LLKA_ConfalPercentile:
            arr[idx].value = float(item)
        else:
            item._fill_c(arr[idx])
    return ctypes.cast(arr, POINTER(c_type)), len(items), arr


class ClassificationContext:
    """Classification parameters needed to classify steps (``LLKA_ClassificationContext``).

    The context may be used as a context manager or closed explicitly with
    :meth:`close` to release the resources early.
    """

    __slots__ = ('_ptr', '_finalizer', '__weakref__')

    def __init__(
        self,
        clusters: Sequence[ClassificationCluster],
        golden_steps: Sequence[GoldenStep],
        confals: Sequence[Confal],
        cluster_nu_angles: Sequence[ClusterNuAngles],
        confal_percentiles: Sequence[float],
        limits: ClassificationLimits,
        max_close_enough_rmsd: float,
    ) -> None:
        """Initializes the context. The context keeps its own copies of all parameters.

        :raises ClassificationContextError: Some of the parameters are invalid. The exception
                                            type identifies the offending parameters.
        """
        c_clusters, n_clusters, _k1 = _to_c_array(clusters, LLKA_ClassificationCluster, 'clusters')
        c_golden, n_golden, _k2 = _to_c_array(golden_steps, LLKA_GoldenStep, 'golden_steps')
        c_confals, n_confals, _k3 = _to_c_array(confals, LLKA_Confal, 'confals')
        c_nus, n_nus, _k4 = _to_c_array(cluster_nu_angles, LLKA_ClusterNuAngles, 'cluster_nu_angles')
        c_pcts, n_pcts, _k5 = _to_c_array(confal_percentiles, LLKA_ConfalPercentile, 'confal_percentiles')
        c_limits = limits._to_c()

        ptr = POINTER(LLKA_ClassificationContext)()
        _capi.check(
            lib.LLKA_initializeClassificationContext(
                c_clusters, n_clusters,
                c_golden, n_golden,
                c_confals, n_confals,
                c_nus, n_nus,
                c_pcts, n_pcts,
                byref(c_limits),
                float(max_close_enough_rmsd),
                byref(ptr),
            ),
            'Cannot initialize classification context',
        )
        self._ptr = ptr
        self._finalizer = weakref.finalize(self, lib.LLKA_destroyClassificationContext, ptr)

    def close(self) -> None:
        """Releases the context. The context cannot be used afterwards."""
        self._finalizer()

    def __enter__(self) -> ClassificationContext:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()

    def _get(self):
        if not self._finalizer.alive:
            raise ValueError('Classification context has been closed')
        return self._ptr

    def classify_step(self, stru: StructureLike) -> ClassifiedStep:
        """Attempts to classify a step.

        Check :attr:`ClassifiedStep.violations` to find out whether the step was fully classified.

        :raises LLKAError: No attempt to classify the structure could be made, e.g. because it is not a valid step.
        """
        stru = as_structure(stru)
        c = LLKA_ClassifiedStep()
        _capi.check(lib.LLKA_classifyStep(byref(stru._c), self._get(), byref(c)), 'Cannot classify step')
        return ClassifiedStep._from_c(c)

    def classify_steps(self, strus: Iterable[StructureLike]) -> list[AttemptedClassifiedStep]:
        """Attempts to classify multiple steps.

        Unlike :meth:`classify_step`, this function does not raise if a step cannot be
        classified; the error is reported in :attr:`AttemptedClassifiedStep.error`.

        :raises NothingToClassifyError: ``strus`` is empty.
        """
        strus = [as_structure(s) for s in strus]
        c_arr = (LLKA_Structure * len(strus))(*[LLKA_Structure(s._c.atoms, s._c.nAtoms) for s in strus])
        c_strus = LLKA_Structures(ctypes.cast(c_arr, POINTER(LLKA_Structure)), len(strus))
        c_results = LLKA_ClassifiedSteps()
        _capi.check(
            lib.LLKA_classifyStepsMultiple(byref(c_strus), self._get(), byref(c_results)),
            'Cannot classify steps',
        )
        try:
            results = []
            for idx in range(c_results.nAttemptedSteps):
                attempt = c_results.attemptedSteps[idx]
                if attempt.status == RetCode.OK:
                    results.append(AttemptedClassifiedStep(ClassifiedStep._from_c(attempt.step), None))
                else:
                    error = error_for_retcode(attempt.status, f'Cannot classify step: {_capi.error_to_string(attempt.status)}')
                    results.append(AttemptedClassifiedStep(None, error))
            return results
        finally:
            lib.LLKA_destroyClassifiedSteps(byref(c_results))

    def average_confal(self, steps: Iterable[ClassifiedStep | AttemptedClassifiedStep]) -> AverageConfal:
        """Calculates the average confal score of classified steps and its percentile.

        Attempted steps that could not be classified are ignored.

        :raises ValueError: The average confal score is not within [0; 100].
        """
        steps = list(steps)
        if all(isinstance(s, AttemptedClassifiedStep) for s in steps) and steps:
            totals = [s.step.confal_score.total for s in steps if s.step is not None]
            self._check_average(totals)

            arr = (LLKA_AttemptedClassifiedStep * len(steps))()
            for idx, s in enumerate(steps):
                if s.step is not None:
                    s.step._fill_c(arr[idx].step)
                    arr[idx].status = RetCode.OK
                else:
                    arr[idx].status = s.error.retcode if s.error is not None else RetCode.E_INVALID_ARGUMENT
            c_steps = LLKA_ClassifiedSteps(ctypes.cast(arr, POINTER(LLKA_AttemptedClassifiedStep)), len(steps))
            c = lib.LLKA_averageConfalAttempted(byref(c_steps), self._get())
        else:
            if not all(isinstance(s, ClassifiedStep) for s in steps):
                raise TypeError('Expected a list of ClassifiedSteps or a list of AttemptedClassifiedSteps')
            self._check_average([s.confal_score.total for s in steps])

            arr = (LLKA_ClassifiedStep * len(steps))()
            for idx, s in enumerate(steps):
                s._fill_c(arr[idx])
            c = lib.LLKA_averageConfal(ctypes.cast(arr, POINTER(LLKA_ClassifiedStep)), len(steps), self._get())

        return AverageConfal(c.score, c.percentile)

    @staticmethod
    def _check_average(totals: list[float]) -> None:
        # libLLKA looks the percentile up without checking the range
        if totals:
            avg = math.fsum(totals) / len(totals)
            if not (0.0 <= avg < 101.0):
                raise ValueError(f'Average confal score {avg} is out of the [0; 100] range')

    def cluster_for_ntc(self, ntc: NtC) -> ClassificationCluster:
        """Returns the classification cluster of a NtC.

        :raises InvalidArgumentError: The context has no cluster for the NtC.
        """
        c = LLKA_ClassificationCluster()
        _capi.check(
            lib.LLKA_classificationClusterForNtC(NtC(ntc), self._get(), byref(c)),
            f'Cannot get classification cluster for {NtC(ntc).name}',
        )
        return ClassificationCluster._from_c(c)

    def confal_for_ntc(self, ntc: NtC) -> Confal:
        """Returns the confal definition of a NtC.

        :raises InvalidArgumentError: The context has no confal for the NtC.
        """
        c = LLKA_Confal()
        _capi.check(
            lib.LLKA_confalForNtC(NtC(ntc), self._get(), byref(c)),
            f'Cannot get confal for {NtC(ntc).name}',
        )
        return Confal._from_c(c)

    def confal_percentile(self, confal_score: float) -> float:
        """Returns the statistical percentile of a confal score.

        :raises ValueError: The score is not within [0; 100].
        """
        confal_score = float(confal_score)
        if not 0.0 <= confal_score <= 100.0:
            raise ValueError(f'Confal score {confal_score} is out of the [0; 100] range')
        return lib.LLKA_confalPercentile(confal_score, self._get())

