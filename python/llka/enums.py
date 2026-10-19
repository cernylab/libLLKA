"""Enumerations mirroring the enums and flag sets of the libLLKA C API."""

from __future__ import annotations

import enum

__all__ = [
    'BaseKind',
    'CANA',
    'ClassificationViolation',
    'CrossResidueMetric',
    'DinucleotideTorsion',
    'NtC',
    'ResidueInStep',
    'ResourceType',
    'RetCode',
    'SugarPucker',
    'SugarPuckerNameBrevity',
    'ViolatingTorsion',
]

_INT32_MAX = 0x7FFFFFFF


class RetCode(enum.IntEnum):
    """Return codes of libLLKA functions (``LLKA_RetCode``)."""

    OK = 0x0
    E_INVALID_ARGUMENT = 0x1
    E_MISMATCHING_SIZES = 0x2
    E_NOT_IMPLEMENTED = 0x3
    E_MULTIPLE_ALT_IDS = 0x4
    E_MISSING_ATOMS = 0x5
    E_MISMATCHING_DATA = 0x6
    E_BAD_CLASSIFICATION_CLUSTERS = 0x7
    E_BAD_GOLDEN_STEPS = 0x8
    E_BAD_CONFALS = 0x9
    E_BAD_AVERAGE_NU_ANGLES = 0x10
    E_BAD_CLASSIFICATION_LIMITS = 0x11
    E_NO_FILE = 0x12
    E_CANNOT_READ_FILE = 0x13
    E_BAD_DATA = 0x14
    E_NO_DATA = 0x15
    E_NOTHING_TO_CLASSIFY = 0x16


class BaseKind(enum.IntEnum):
    """Kind of a nucleobase (``LLKA_BaseKind``)."""

    PURINE = 0
    PYRIMIDINE = 1
    NON_STANDARD_BASE = 2


class CrossResidueMetric(enum.IntEnum):
    """NtC step metrics that span both residues of a step (``LLKA_CrossResidueMetric``)."""

    DIST_CC = 0
    DIST_NN = 1
    TOR_MU = 2


class DinucleotideTorsion(enum.IntEnum):
    """Backbone and base torsions of a dinucleotide (``LLKA_DinucleotideTorsion``)."""

    DELTA_1 = 0
    EPSILON_1 = 1
    ZETA_1 = 2
    ALPHA_2 = 3
    BETA_2 = 4
    GAMMA_2 = 5
    DELTA_2 = 6
    CHI_1 = 7
    CHI_2 = 8


class ResidueInStep(enum.IntEnum):
    """Position of a residue in a dinucleotide step (``LLKA_ResidueInStep``)."""

    FIRST_RESIDUE = 0
    SECOND_RESIDUE = 1


class CANA(enum.IntEnum):
    """CANA classes (``LLKA_CANA``).

    ``INVALID`` is used by the C API to indicate "no CANA". The Python API
    reports a missing CANA as ``None`` instead.
    """

    AAA = 0
    AAw = 1
    AAu = 2
    AB = 3
    BA = 4
    BBB = 5
    BBw = 6
    B12 = 7
    BB2 = 8
    miB = 9
    ICL = 10
    OPN = 11
    SYN = 12
    ZZZ = 13
    INVALID = _INT32_MAX


class NtC(enum.IntEnum):
    """NtC classes (``LLKA_NtC``).

    ``INVALID`` is used by the C API to indicate "no NtC". The Python API
    reports a missing NtC as ``None`` instead.
    """

    AA00 = 0
    AA02 = 1
    AA03 = 2
    AA04 = 3
    AA08 = 4
    AA09 = 5
    AA01 = 6
    AA05 = 7
    AA06 = 8
    AA10 = 9
    AA11 = 10
    AA07 = 11
    AA12 = 12
    AA13 = 13
    AB01 = 14
    AB02 = 15
    AB03 = 16
    AB04 = 17
    AB05 = 18
    BA01 = 19
    BA05 = 20
    BA09 = 21
    BA08 = 22
    BA10 = 23
    BA13 = 24
    BA16 = 25
    BA17 = 26
    BB00 = 27
    BB01 = 28
    BB17 = 29
    BB02 = 30
    BB03 = 31
    BB11 = 32
    BB16 = 33
    BB04 = 34
    BB05 = 35
    BB07 = 36
    BB08 = 37
    BB10 = 38
    BB12 = 39
    BB13 = 40
    BB14 = 41
    BB15 = 42
    BB20 = 43
    IC01 = 44
    IC02 = 45
    IC03 = 46
    IC04 = 47
    IC05 = 48
    IC06 = 49
    IC07 = 50
    OP01 = 51
    OP02 = 52
    OP03 = 53
    OP04 = 54
    OP05 = 55
    OP06 = 56
    OP07 = 57
    OP08 = 58
    OP09 = 59
    OP10 = 60
    OP11 = 61
    OP12 = 62
    OP13 = 63
    OP14 = 64
    OP15 = 65
    OP16 = 66
    OP17 = 67
    OP18 = 68
    OP19 = 69
    OP20 = 70
    OP21 = 71
    OP22 = 72
    OP23 = 73
    OP24 = 74
    OP25 = 75
    OP26 = 76
    OP27 = 77
    OP28 = 78
    OP29 = 79
    OP30 = 80
    OP31 = 81
    OPS1 = 82
    OP1S = 83
    AAS1 = 84
    AB1S = 85
    AB2S = 86
    BB1S = 87
    BB2S = 88
    BBS1 = 89
    ZZ01 = 90
    ZZ02 = 91
    ZZ1S = 92
    ZZ2S = 93
    ZZS1 = 94
    ZZS2 = 95
    INVALID = _INT32_MAX


class SugarPucker(enum.IntEnum):
    """Sugar pucker conformations (``LLKA_SugarPucker``)."""

    C3_ENDO = 0
    C4_EXO = 1
    O4_ENDO = 2
    C1_EXO = 3
    C2_ENDO = 4
    C3_EXO = 5
    C4_ENDO = 6
    O4_EXO = 7
    C1_ENDO = 8
    C2_EXO = 9
    INVALID = _INT32_MAX


class SugarPuckerNameBrevity(enum.IntEnum):
    """Level of brevity of a sugar pucker name (``LLKA_SugarPuckerNameBrevity``)."""

    VERY_TERSE = 0  # C1exo, C1end
    TERSE = 1       # C1exo, C1endo
    FANCY = 2       # "C1' exo", "C1' endo"


class ClassificationViolation(enum.IntFlag):
    """Reasons that prevented assignment of NtC and CANA classes to a step."""

    OK = 0
    SCORE_TOO_LOW = 1 << 0
    NOT_ENOUGH_NEAREST_NEIGHBORS = 1 << 1
    AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT = 1 << 2
    NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT = 1 << 3
    CC_TOO_LOW = 1 << 4
    CC_TOO_HIGH = 1 << 5
    NN_TOO_LOW = 1 << 6
    NN_TOO_HIGH = 1 << 7
    MU_TOO_LOW = 1 << 8
    MU_TOO_HIGH = 1 << 9
    TOTAL_DISTANCE_TOO_HIGH = 1 << 10
    FIRST_PSEUDOROTATION_TOO_DIFFERENT = 1 << 11
    SECOND_PSEUDOROTATION_TOO_DIFFERENT = 1 << 12
    BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES = 1 << 13
    DELTA_TORSION_ANGLE_REJECTED = 1 << 14
    WRONG_METRICS = 1 << 15
    UNASSIGNED_BUT_CLOSE_ENOUGH = 1 << 16


class ViolatingTorsion(enum.IntFlag):
    """Torsions that exceeded classification tolerances."""

    NONE = 0
    DELTA_1 = 1 << 0
    EPSILON_1 = 1 << 1
    ZETA_1 = 1 << 2
    ALPHA_2 = 1 << 3
    BETA_2 = 1 << 4
    GAMMA_2 = 1 << 5
    DELTA_2 = 1 << 6
    CHI_1 = 1 << 7
    CHI_2 = 1 << 8


class ResourceType(enum.IntEnum):
    """Kinds of classification parameter files that libLLKA can load (``LLKA_ResourceType``)."""

    AVERAGE_NU_ANGLES = 0
    CLUSTERS = 1
    CONFALS = 2
    GOLDEN_STEPS = 3
    CONFAL_PERCENTILES = 4
