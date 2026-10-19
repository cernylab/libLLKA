"""Low-level ctypes interface to the libLLKA C API.

This module loads the libLLKA shared library, mirrors the C data structures
as ctypes structures and declares the prototypes of the C functions. It is an
implementation detail of the ``llka`` package; use the public modules instead.

The library is looked up in this order:

1. The path given by the ``LLKA_LIBRARY`` environment variable.
2. The directory of this package (installed packages ship the library here).
3. The standard library search path of the platform.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import struct
import sys
from ctypes import (
    POINTER,
    Structure,
    c_char,
    c_char_p,
    c_double,
    c_float,
    c_int16,
    c_int32,
    c_size_t,
    c_ubyte,
    c_uint32,
    c_void_p,
    c_wchar_p,
)
from pathlib import Path

from .errors import error_for_retcode

# Enums are passed as int32. Enums that are not forced to int32 in the C headers
# only have small values and are therefore int-sized as well.
c_enum = c_int32
LLKA_Bool = c_ubyte
LLKA_RetCode = c_enum

LLKA_PathChar_p = c_wchar_p if sys.platform == 'win32' else c_char_p

LLKA_NO_ALTID = b'\0'
LLKA_INVALID_BKBN_ATOM_INDEX = ctypes.c_size_t(-1).value

#
# llka_main.h
#

class LLKA_Point(Structure):
    _fields_ = [('x', c_double), ('y', c_double), ('z', c_double)]


class LLKA_Points(Structure):
    # The C struct wraps the pointer in a union of LLKA_Point * and double *
    _fields_ = [('points', POINTER(LLKA_Point)), ('nPoints', c_size_t)]


class LLKA_Matrix(Structure):
    _fields_ = [('data', POINTER(c_double)), ('nCols', c_size_t), ('nRows', c_size_t)]


#
# llka_structure.h
#

class LLKA_AlternatePositions(Structure):
    _fields_ = [('positions', POINTER(c_char)), ('nPositions', c_size_t)]


class LLKA_Atom(Structure):
    _fields_ = [
        ('type_symbol', c_char_p),
        ('label_atom_id', c_char_p),
        ('label_entity_id', c_char_p),
        ('label_comp_id', c_char_p),
        ('label_asym_id', c_char_p),
        ('auth_atom_id', c_char_p),
        ('auth_comp_id', c_char_p),
        ('auth_asym_id', c_char_p),
        ('coords', LLKA_Point),
        ('id', c_uint32),
        ('label_seq_id', c_int32),
        ('auth_seq_id', c_int32),
        ('pdbx_PDB_model_num', c_int32),
        ('pdbx_PDB_ins_code', c_char_p),
        ('label_alt_id', c_char),
    ]


class LLKA_Structure(Structure):
    _fields_ = [('atoms', POINTER(LLKA_Atom)), ('nAtoms', c_size_t)]


class LLKA_StructureView(Structure):
    _fields_ = [('atoms', POINTER(POINTER(LLKA_Atom))), ('nAtoms', c_size_t), ('capacity', c_size_t)]


class LLKA_Structures(Structure):
    _fields_ = [('strus', POINTER(LLKA_Structure)), ('nStrus', c_size_t)]


#
# llka_ntc.h
#

class LLKA_AtomNameQuad(Structure):
    _fields_ = [('a', c_char * 8), ('b', c_char * 8), ('c', c_char * 8), ('d', c_char * 8)]


class LLKA_BackboneAtom(Structure):
    _fields_ = [('residue', c_enum), ('name', c_char_p)]


class LLKA_StepInfo(Structure):
    _fields_ = [
        ('firstSeqId', c_int32),
        ('firstSeqIdAuth', c_int32),
        ('secondSeqId', c_int32),
        ('secondSeqIdAuth', c_int32),
        ('firstBaseKind', c_enum),
        ('secondBaseKind', c_enum),
    ]


STEP_METRICS_FIELDS = (
    'delta_1', 'epsilon_1', 'zeta_1', 'alpha_2', 'beta_2', 'gamma_2', 'delta_2', 'chi_1', 'chi_2', 'CC', 'NN', 'mu'
)


class LLKA_StepMetrics(Structure):
    _fields_ = [(name, c_double) for name in STEP_METRICS_FIELDS]


#
# llka_nucleotide.h
#

NU_ANGLES_FIELDS = ('nu_0', 'nu_1', 'nu_2', 'nu_3', 'nu_4')


class LLKA_NuAngles(Structure):
    _fields_ = [(name, c_double) for name in NU_ANGLES_FIELDS]


class LLKA_RiboseMetrics(Structure):
    _fields_ = [('nus', LLKA_NuAngles), ('P', c_double), ('tMax', c_double), ('pucker', c_enum)]


#
# llka_segmentation.h
#

class LLKA_Residue(Structure):
    _fields_ = [('atoms', POINTER(POINTER(LLKA_Atom))), ('nAtoms', c_size_t)]


class LLKA_Chain(Structure):
    _fields_ = [
        ('residues', POINTER(LLKA_Residue)),
        ('nResidues', c_size_t),
        ('atoms', POINTER(POINTER(LLKA_Atom))),
        ('nAtoms', c_size_t),
    ]


class LLKA_Model(Structure):
    _fields_ = [
        ('chains', POINTER(LLKA_Chain)),
        ('nChains', c_size_t),
        ('atoms', POINTER(POINTER(LLKA_Atom))),
        ('nAtoms', c_size_t),
    ]


class LLKA_StructureSegments(Structure):
    _fields_ = [('models', POINTER(LLKA_Model)), ('nModels', c_size_t)]


#
# llka_connectivity_similarity.h
#

class LLKA_Connectivity(Structure):
    _fields_ = [('C5PrimeDistance', c_double), ('O3PrimeDistance', c_double)]


class LLKA_Connectivities(Structure):
    _fields_ = [('conns', POINTER(LLKA_Connectivity)), ('nConns', c_size_t)]


class LLKA_Similarity(Structure):
    _fields_ = [('rmsd', c_double), ('euclideanDistance', c_double)]


class LLKA_Similarities(Structure):
    _fields_ = [('similars', POINTER(LLKA_Similarity)), ('nSimilars', c_size_t)]


#
# llka_classification.h
#

class LLKA_AverageConfal(Structure):
    _fields_ = [('score', c_double), ('percentile', c_double)]


class LLKA_ClassificationLimits(Structure):
    _fields_ = [
        ('minimumNearestNeighbors', c_size_t),
        ('numberOfUsedNearestNeighbors', c_size_t),
        ('averageNeighborsTorsionCutoff', c_double),
        ('nearestNeighborTorsionsCutoff', c_double),
        ('totalDistanceCutoff', c_double),
        ('pseudorotationCutoff', c_double),
        ('minimumClusterVotes', c_double),
    ]


class LLKA_ClassificationMetric(Structure):
    _fields_ = [('deviation', c_double), ('minValue', c_double), ('meanValue', c_double), ('maxValue', c_double)]


class LLKA_NuAnglesMetrics(Structure):
    _fields_ = [(name, LLKA_ClassificationMetric) for name in NU_ANGLES_FIELDS]


class LLKA_ClusterNuAngles(Structure):
    _fields_ = [
        ('firstNucleotide', LLKA_NuAnglesMetrics),
        ('secondNucleotide', LLKA_NuAnglesMetrics),
        ('clusterNumber', c_int32),
    ]


class LLKA_ClassificationCluster(Structure):
    _fields_ = [(name, LLKA_ClassificationMetric) for name in STEP_METRICS_FIELDS] + [
        ('nusFirst', LLKA_NuAnglesMetrics),
        ('nusSecond', LLKA_NuAnglesMetrics),
        ('ribosePseudorotation_1', c_double),
        ('ribosePseudorotation_2', c_double),
        ('number', c_int32),
        ('NtC', c_enum),
        ('CANA', c_enum),
    ]


class LLKA_ConfalPercentile(Structure):
    _fields_ = [('value', c_double)]


class LLKA_ConfalScore(Structure):
    _fields_ = [(name, c_double) for name in STEP_METRICS_FIELDS] + [('total', c_double)]


class LLKA_ClassifiedStep(Structure):
    _fields_ = [
        ('assignedNtC', c_enum),
        ('assignedCANA', c_enum),
        ('closestNtC', c_enum),
        ('closestCANA', c_enum),
        ('confalScore', LLKA_ConfalScore),
        ('euclideanDistanceNtCIdeal', c_double),
        ('metrics', LLKA_StepMetrics),
        ('differencesFromNtCAverages', LLKA_StepMetrics),
        ('nuAngles_1', LLKA_NuAngles),
        ('nuAngles_2', LLKA_NuAngles),
        ('ribosePseudorotation_1', c_double),
        ('ribosePseudorotation_2', c_double),
        ('tau_1', c_double),
        ('tau_2', c_double),
        ('sugarPucker_1', c_enum),
        ('sugarPucker_2', c_enum),
        ('nuAngleDifferences_1', LLKA_NuAngles),
        ('nuAngleDifferences_2', LLKA_NuAngles),
        ('rmsdToClosestNtC', c_double),
        ('closestGoldenStep', c_char_p),
        ('violations', c_int32),
        ('violatingTorsionsAverage', c_int16),
        ('violatingTorsionsNearest', c_int16),
    ]


class LLKA_AttemptedClassifiedStep(Structure):
    _fields_ = [('step', LLKA_ClassifiedStep), ('status', LLKA_RetCode)]


class LLKA_ClassifiedSteps(Structure):
    _fields_ = [('attemptedSteps', POINTER(LLKA_AttemptedClassifiedStep)), ('nAttemptedSteps', c_size_t)]


class LLKA_Confal(Structure):
    _fields_ = [(name, c_double) for name in STEP_METRICS_FIELDS] + [
        ('nusFirst', LLKA_NuAngles),
        ('nusSecond', LLKA_NuAngles),
        ('clusterNumber', c_int32),
    ]


class LLKA_GoldenStep(Structure):
    _fields_ = [
        ('pucker_1', c_enum),
        ('pucker_2', c_enum),
        ('nuAngles_1', LLKA_NuAngles),
        ('nuAngles_2', LLKA_NuAngles),
        ('metrics', LLKA_StepMetrics),
        ('name', c_char_p),
        ('clusterIdx', c_size_t),
        ('clusterNumber', c_int32),
    ]


class LLKA_ClassificationContext(Structure):
    """Opaque type, only ever handled through a pointer."""


#
# llka_resource_loaders.h
#

class LLKA_Resource(Structure):
    # "data" is a union of pointers to the arrays of the individual resource types
    _fields_ = [('type', c_enum), ('data', c_void_p), ('count', c_size_t)]


#
# llka_minicif.h
#

class LLKA_CifDataValue(Structure):
    _fields_ = [('text', c_char_p), ('state', c_enum)]


class LLKA_CifDataItem(Structure):
    pass


LLKA_CifDataItem._fields_ = [
    ('keyword', c_char_p),
    ('values', POINTER(LLKA_CifDataValue)),
    ('nValues', c_size_t),
    ('p', c_void_p),
]


class LLKA_CifDataCategory(Structure):
    _fields_ = [('name', c_char_p), ('firstItem', POINTER(LLKA_CifDataItem)), ('p', c_void_p)]


class LLKA_CifDataBlock(Structure):
    _fields_ = [('name', c_char_p), ('firstCategory', POINTER(LLKA_CifDataCategory)), ('p', c_void_p)]


class LLKA_CifData(Structure):
    _fields_ = [('blocks', POINTER(LLKA_CifDataBlock)), ('nBlocks', c_size_t), ('p', c_void_p)]


class LLKA_StructureEntry(Structure):
    _fields_ = [('id', c_char_p)]


class LLKA_ImportedStructure(Structure):
    _fields_ = [('entry', LLKA_StructureEntry), ('structure', LLKA_Structure), ('cifData', POINTER(LLKA_CifData))]


#
# llka_tracing.h
#

class LLKA_TracepointInfo(Structure):
    _fields_ = [('TPID', c_int32), ('description', c_char_p)]


class LLKA_TracepointInfos(Structure):
    _fields_ = [('infos', POINTER(LLKA_TracepointInfo)), ('nInfos', c_size_t)]


# All ctypes mirrors of C structures. Used to verify that the layouts match the C headers.
STRUCTURES = {
    cls.__name__: cls
    for cls in (
        LLKA_Point, LLKA_Points, LLKA_Matrix,
        LLKA_AlternatePositions, LLKA_Atom, LLKA_Structure, LLKA_StructureView, LLKA_Structures,
        LLKA_AtomNameQuad, LLKA_BackboneAtom, LLKA_StepInfo, LLKA_StepMetrics,
        LLKA_NuAngles, LLKA_RiboseMetrics,
        LLKA_Residue, LLKA_Chain, LLKA_Model, LLKA_StructureSegments,
        LLKA_Connectivity, LLKA_Connectivities, LLKA_Similarity, LLKA_Similarities,
        LLKA_AverageConfal, LLKA_ClassificationLimits, LLKA_ClassificationMetric, LLKA_NuAnglesMetrics,
        LLKA_ClusterNuAngles, LLKA_ClassificationCluster, LLKA_ConfalPercentile, LLKA_ConfalScore,
        LLKA_ClassifiedStep, LLKA_AttemptedClassifiedStep, LLKA_ClassifiedSteps, LLKA_Confal, LLKA_GoldenStep,
        LLKA_Resource,
        LLKA_CifDataValue, LLKA_CifDataItem, LLKA_CifDataCategory, LLKA_CifDataBlock, LLKA_CifData,
        LLKA_StructureEntry, LLKA_ImportedStructure,
        LLKA_TracepointInfo, LLKA_TracepointInfos,
    )
}

_P = POINTER

# Function prototypes: name -> (restype, argtypes).
#
# Functions that are declared in the public headers but cannot be called through
# the C ABI are listed here and left out of the prototypes.
BROKEN_FUNCTIONS: frozenset[str] = frozenset()

# Functions whose prototypes are declared but that the Python API does not need.
UNWRAPPED_FUNCTIONS = {
    'LLKA_appendAtomFromParams': 'Structure.append() takes an Atom value',
    'LLKA_duplicateAtom': 'Atom is an immutable Python value',
    'LLKA_emptyAtom': 'Atom is an immutable Python value',
    'LLKA_makeAtom': 'Atom is an immutable Python value',
    'LLKA_cifDataBlock_initialize': 'only needed to initialize blocks preallocated by C code',
}

# Functions that are not available when libLLKA is built with LLKA_FILESYSTEM_ACCESS_DISABLED
FILESYSTEM_FUNCTIONS = frozenset({
    'LLKA_cifFileToData',
    'LLKA_cifFileToStructure',
    'LLKA_loadResourceFile',
})

PROTOTYPES: dict[str, tuple[object, list[object]]] = {
    # llka_main.h
    'LLKA_destroyMatrix': (None, [_P(LLKA_Matrix)]),
    'LLKA_destroyPoints': (None, [_P(LLKA_Points)]),
    'LLKA_destroyString': (None, [c_void_p]),
    'LLKA_errorToString': (c_char_p, [LLKA_RetCode]),
    'LLKA_initMatrix': (None, [c_size_t, c_size_t, _P(LLKA_Matrix)]),
    'LLKA_duplicateMatrix': (LLKA_Matrix, [_P(LLKA_Matrix)]),
    'LLKA_makePoints': (LLKA_Points, [_P(LLKA_Point), c_size_t]),
    'LLKA_zeroPoints': (LLKA_Points, [c_size_t]),

    # llka_structure.h
    'LLKA_appendAtom': (None, [_P(LLKA_Atom), _P(LLKA_Structure)]),
    'LLKA_appendAtomFromParams': (None, [
        c_uint32, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p,
        c_int32, c_char, c_int32, c_char_p, c_int32, _P(LLKA_Point), _P(LLKA_Structure),
    ]),
    'LLKA_appendStructure': (None, [_P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_baseKind': (LLKA_RetCode, [c_char_p, _P(c_enum)]),
    'LLKA_compareAtoms': (LLKA_Bool, [_P(LLKA_Atom), _P(LLKA_Atom), LLKA_Bool]),
    'LLKA_destroyAlternatePositions': (None, [_P(LLKA_AlternatePositions)]),
    'LLKA_destroyAtom': (None, [_P(LLKA_Atom)]),
    'LLKA_destroyStructure': (None, [_P(LLKA_Structure)]),
    'LLKA_destroyStructureView': (None, [_P(LLKA_StructureView)]),
    'LLKA_destroyStructures': (None, [_P(LLKA_Structures)]),
    'LLKA_duplicateAtom': (None, [_P(LLKA_Atom), _P(LLKA_Atom)]),
    'LLKA_duplicateStructure': (None, [_P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_emptyAtom': (LLKA_Atom, []),
    'LLKA_findAtom': (_P(LLKA_Atom), [
        _P(LLKA_Structure), c_char_p, c_char_p, c_char_p, c_int32, c_char, c_char_p, c_int32,
    ]),
    'LLKA_findAtomView': (_P(LLKA_Atom), [
        _P(LLKA_StructureView), c_char_p, c_char_p, c_char_p, c_int32, c_char, c_char_p, c_int32,
    ]),
    'LLKA_initAtom': (None, [
        c_uint32, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p,
        c_int32, c_char, c_int32, c_char_p, c_int32, _P(LLKA_Point), _P(LLKA_Atom),
    ]),
    'LLKA_initStructure': (None, [_P(LLKA_Atom), c_size_t, _P(LLKA_Structure)]),
    'LLKA_makeAtom': (LLKA_Atom, [
        c_uint32, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p,
        c_int32, c_char, c_int32, c_char_p, c_int32, _P(LLKA_Point),
    ]),
    'LLKA_makeStructure': (LLKA_Structure, [_P(LLKA_Atom), c_size_t]),
    'LLKA_makeStructureFromPtrs': (LLKA_Structure, [_P(_P(LLKA_Atom)), c_size_t]),
    'LLKA_removeAtomById': (None, [c_uint32, _P(LLKA_Structure)]),
    'LLKA_splitByAltIds': (LLKA_Structures, [_P(LLKA_Structure), _P(LLKA_AlternatePositions)]),
    'LLKA_splitStructureToDinucleotideSteps': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structures)]),

    # llka_ntc.h
    'LLKA_backboneAtomIndex': (c_size_t, [_P(LLKA_BackboneAtom), _P(LLKA_Structure)]),
    'LLKA_CANAToName': (c_char_p, [c_enum]),
    'LLKA_calculateStepMetrics': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_StepMetrics)]),
    'LLKA_calculateStepMetricsDifferenceAgainstReference': (
        LLKA_RetCode, [_P(LLKA_Structure), c_enum, _P(LLKA_StepMetrics)]
    ),
    'LLKA_crossResidueMetric': (LLKA_RetCode, [c_enum, _P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_crossResidueMetricAtomsFromBases': (LLKA_RetCode, [c_char_p, c_char_p, c_enum, _P(LLKA_AtomNameQuad)]),
    'LLKA_crossResidueMetricAtomsFromStructure': (LLKA_RetCode, [_P(LLKA_Structure), c_enum, _P(LLKA_AtomNameQuad)]),
    'LLKA_crossResidueMetricName': (c_char_p, [c_enum, LLKA_Bool]),
    'LLKA_dinucleotideTorsion': (LLKA_RetCode, [c_enum, _P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_dinucleotideTorsionAtomsFromBases': (LLKA_RetCode, [c_char_p, c_char_p, c_enum, _P(LLKA_AtomNameQuad)]),
    'LLKA_dinucleotideTorsionAtomsFromStructure': (
        LLKA_RetCode, [_P(LLKA_Structure), c_enum, _P(LLKA_AtomNameQuad)]
    ),
    'LLKA_dinucleotideTorsionName': (c_char_p, [c_enum, LLKA_Bool]),
    'LLKA_extractBackbone': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_extractBackboneView': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_StructureView)]),
    'LLKA_extractExtendedBackbone': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_extractExtendedBackboneView': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_StructureView)]),
    'LLKA_extractMetricsStructure': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_nameToCANA': (c_enum, [c_char_p]),
    'LLKA_nameToNtC': (c_enum, [c_char_p]),
    'LLKA_NtCStructure': (LLKA_Structure, [c_enum]),
    'LLKA_NtCToName': (c_char_p, [c_enum]),
    'LLKA_structureIsStep': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_StepInfo)]),

    # llka_measurements.h
    'LLKA_measureAnglef': (c_float, [_P(LLKA_Atom), _P(LLKA_Atom), _P(LLKA_Atom)]),
    'LLKA_measureAngle': (c_double, [_P(LLKA_Atom), _P(LLKA_Atom), _P(LLKA_Atom)]),
    'LLKA_measureDihedralf': (c_float, [_P(LLKA_Atom), _P(LLKA_Atom), _P(LLKA_Atom), _P(LLKA_Atom)]),
    'LLKA_measureDihedral': (c_double, [_P(LLKA_Atom), _P(LLKA_Atom), _P(LLKA_Atom), _P(LLKA_Atom)]),
    'LLKA_measureDistancef': (c_float, [_P(LLKA_Atom), _P(LLKA_Atom)]),
    'LLKA_measureDistance': (c_double, [_P(LLKA_Atom), _P(LLKA_Atom)]),

    # llka_nucleotide.h
    'LLKA_extractNucleotide': (LLKA_Structure, [_P(LLKA_Structure), c_int32, c_char_p, c_int32]),
    'LLKA_extractNucleotideView': (LLKA_StructureView, [_P(LLKA_Structure), c_int32, c_char_p, c_int32]),
    'LLKA_extractRibose': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure)]),
    'LLKA_extractRiboseView': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_StructureView)]),
    'LLKA_isNucleotideCompound': (LLKA_Bool, [c_char_p]),
    'LLKA_nameToSugarPucker': (c_enum, [c_char_p]),
    'LLKA_riboseMetrics': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_RiboseMetrics)]),
    'LLKA_riboseMetricsView': (LLKA_RetCode, [_P(LLKA_StructureView), _P(LLKA_RiboseMetrics)]),
    'LLKA_sugarPucker': (LLKA_RetCode, [_P(LLKA_Structure), _P(c_enum)]),
    'LLKA_sugarPuckerView': (LLKA_RetCode, [_P(LLKA_StructureView), _P(c_enum)]),
    'LLKA_sugarPuckerToName': (c_char_p, [c_enum, c_enum]),

    # llka_superposition.h
    'LLKA_applyTransformationPoints': (LLKA_RetCode, [_P(LLKA_Points), _P(LLKA_Matrix)]),
    'LLKA_applyTransformationStructure': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Matrix)]),
    'LLKA_centroidPoints': (LLKA_Point, [_P(LLKA_Points)]),
    'LLKA_centroidStructure': (LLKA_Point, [_P(LLKA_Structure)]),
    'LLKA_rmsdPoints': (LLKA_RetCode, [_P(LLKA_Points), _P(LLKA_Points), _P(c_double)]),
    'LLKA_rmsdStructures': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure), _P(c_double)]),
    'LLKA_superposePoints': (LLKA_RetCode, [_P(LLKA_Points), _P(LLKA_Points), _P(c_double)]),
    'LLKA_superposeStructures': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure), _P(c_double)]),
    'LLKA_superposeStructuresView': (LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_StructureView), _P(c_double)]),
    'LLKA_superpositionMatrixPoints': (LLKA_RetCode, [_P(LLKA_Points), _P(LLKA_Points), _P(LLKA_Matrix)]),
    'LLKA_superpositionMatrixStructures': (
        LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Matrix)]
    ),
    'LLKA_superpositionMatrixStructureViews': (
        LLKA_RetCode, [_P(LLKA_StructureView), _P(LLKA_StructureView), _P(LLKA_Matrix)]
    ),

    # llka_segmentation.h
    'LLKA_destroyStructureSegments': (None, [_P(LLKA_StructureSegments)]),
    'LLKA_structureSegments': (LLKA_StructureSegments, [_P(LLKA_Structure)]),

    # llka_connectivity_similarity.h
    'LLKA_measureStepConnectivityNtCs': (
        LLKA_RetCode, [_P(LLKA_Structure), c_enum, _P(LLKA_Structure), c_enum, _P(LLKA_Connectivity)]
    ),
    'LLKA_measureStepConnectivityNtCsMultipleFirst': (
        LLKA_RetCode, [_P(LLKA_Structure), _P(c_enum), _P(LLKA_Structure), c_enum, _P(LLKA_Connectivities)]
    ),
    'LLKA_measureStepConnectivityNtCsMultipleSecond': (
        LLKA_RetCode, [_P(LLKA_Structure), c_enum, _P(LLKA_Structure), _P(c_enum), _P(LLKA_Connectivities)]
    ),
    'LLKA_measureStepConnectivityStructures': (LLKA_RetCode, [
        _P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Connectivity),
    ]),
    'LLKA_measureStepConnectivityStructuresMultiple': (LLKA_RetCode, [
        _P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Structures), _P(LLKA_Connectivities),
    ]),
    'LLKA_measureStepSimilarityNtC': (LLKA_RetCode, [_P(LLKA_Structure), c_enum, _P(LLKA_Similarity)]),
    'LLKA_measureStepSimilarityNtCMultiple': (LLKA_RetCode, [_P(LLKA_Structure), _P(c_enum), _P(LLKA_Similarities)]),
    'LLKA_measureStepSimilarityStructure': (
        LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_Structure), _P(LLKA_Similarity)]
    ),

    # llka_classification.h
    'LLKA_averageConfal': (LLKA_AverageConfal, [_P(LLKA_ClassifiedStep), c_size_t, _P(LLKA_ClassificationContext)]),
    'LLKA_averageConfalAttempted': (
        LLKA_AverageConfal, [_P(LLKA_ClassifiedSteps), _P(LLKA_ClassificationContext)]
    ),
    'LLKA_classificationClusterForNtC': (
        LLKA_RetCode, [c_enum, _P(LLKA_ClassificationContext), _P(LLKA_ClassificationCluster)]
    ),
    'LLKA_classificationViolationToName': (c_char_p, [c_int32]),
    'LLKA_classifyStep': (
        LLKA_RetCode, [_P(LLKA_Structure), _P(LLKA_ClassificationContext), _P(LLKA_ClassifiedStep)]
    ),
    'LLKA_classifyStepsMultiple': (
        LLKA_RetCode, [_P(LLKA_Structures), _P(LLKA_ClassificationContext), _P(LLKA_ClassifiedSteps)]
    ),
    'LLKA_confalForNtC': (LLKA_RetCode, [c_enum, _P(LLKA_ClassificationContext), _P(LLKA_Confal)]),
    'LLKA_confalPercentile': (c_double, [c_double, _P(LLKA_ClassificationContext)]),
    'LLKA_destroyClassificationContext': (None, [_P(LLKA_ClassificationContext)]),
    'LLKA_destroyClassifiedSteps': (None, [_P(LLKA_ClassifiedSteps)]),
    'LLKA_initializeClassificationContext': (LLKA_RetCode, [
        _P(LLKA_ClassificationCluster), c_size_t,
        _P(LLKA_GoldenStep), c_size_t,
        _P(LLKA_Confal), c_size_t,
        _P(LLKA_ClusterNuAngles), c_size_t,
        _P(LLKA_ConfalPercentile), c_size_t,
        _P(LLKA_ClassificationLimits),
        c_double,
        _P(_P(LLKA_ClassificationContext)),
    ]),

    # llka_resource_loaders.h
    'LLKA_destroyResource': (None, [_P(LLKA_Resource)]),
    'LLKA_loadResourceFile': (LLKA_RetCode, [LLKA_PathChar_p, _P(LLKA_Resource)]),
    'LLKA_loadResourceText': (LLKA_RetCode, [c_char_p, _P(LLKA_Resource)]),

    # llka_minicif.h
    'LLKA_cifData_addBlock': (None, [_P(LLKA_CifData), c_char_p]),
    'LLKA_cifData_detaint': (LLKA_RetCode, [_P(LLKA_CifData)]),
    'LLKA_cifData_deleteBlock': (LLKA_RetCode, [_P(LLKA_CifData), c_size_t]),
    'LLKA_cifData_empty': (_P(LLKA_CifData), []),
    'LLKA_cifData_isTainted': (LLKA_Bool, [_P(LLKA_CifData)]),
    'LLKA_cifDataBlock_addCategory': (_P(LLKA_CifDataCategory), [_P(LLKA_CifDataBlock), c_char_p]),
    'LLKA_cifDataBlock_deleteCategory': (LLKA_RetCode, [_P(LLKA_CifDataBlock), c_char_p]),
    'LLKA_cifDataBlock_findCategory': (_P(LLKA_CifDataCategory), [_P(LLKA_CifDataBlock), c_char_p]),
    'LLKA_cifDataBlock_initialize': (LLKA_RetCode, [_P(LLKA_CifDataBlock), c_char_p, _P(LLKA_CifData)]),
    'LLKA_cifDataBlock_nextCategory': (_P(LLKA_CifDataCategory), [_P(LLKA_CifDataCategory)]),
    'LLKA_cifDataBlock_previousCategory': (_P(LLKA_CifDataCategory), [_P(LLKA_CifDataCategory)]),
    'LLKA_cifDataCategory_addItem': (_P(LLKA_CifDataItem), [_P(LLKA_CifDataCategory), c_char_p]),
    'LLKA_cifDataCategory_deleteItem': (LLKA_RetCode, [_P(LLKA_CifDataCategory), c_char_p]),
    'LLKA_cifDataCategory_findItem': (_P(LLKA_CifDataItem), [_P(LLKA_CifDataCategory), c_char_p]),
    'LLKA_cifDataCategory_isLoop': (LLKA_Bool, [_P(LLKA_CifDataCategory)]),
    'LLKA_cifDataCategory_nextItem': (_P(LLKA_CifDataItem), [_P(LLKA_CifDataItem)]),
    'LLKA_cifDataCategory_previousItem': (_P(LLKA_CifDataItem), [_P(LLKA_CifDataItem)]),
    'LLKA_cifDataItem_addValue': (None, [_P(LLKA_CifDataItem), _P(LLKA_CifDataValue)]),
    'LLKA_cifDataItem_addValues': (None, [_P(LLKA_CifDataItem), _P(LLKA_CifDataValue), c_size_t]),
    'LLKA_cifDataItem_setValues': (None, [_P(LLKA_CifDataItem), _P(LLKA_CifDataValue), c_size_t]),
    'LLKA_cifDataToString': (LLKA_RetCode, [_P(LLKA_CifData), LLKA_Bool, _P(c_void_p)]),
    'LLKA_cifFileToData': (LLKA_RetCode, [LLKA_PathChar_p, _P(_P(LLKA_CifData)), _P(c_void_p)]),
    'LLKA_cifFileToStructure': (
        LLKA_RetCode, [LLKA_PathChar_p, _P(LLKA_ImportedStructure), _P(c_void_p), c_int32]
    ),
    'LLKA_cifTextToData': (LLKA_RetCode, [c_char_p, _P(_P(LLKA_CifData)), _P(c_void_p)]),
    'LLKA_cifTextToStructure': (LLKA_RetCode, [c_char_p, _P(LLKA_ImportedStructure), _P(c_void_p), c_int32]),
    'LLKA_destroyCifData': (None, [_P(LLKA_CifData)]),
    'LLKA_destroyImportedStructure': (None, [_P(LLKA_ImportedStructure)]),

    # llka_tracing.h
    'LLKA_destroyTrace': (None, [c_void_p]),
    'LLKA_destroyTracepointInfo': (None, [_P(LLKA_TracepointInfos)]),
    'LLKA_toggleAllTracepoints': (None, [LLKA_Bool]),
    'LLKA_toggleTracepoint': (None, [c_int32, LLKA_Bool]),
    'LLKA_trace': (c_void_p, [LLKA_Bool]),
    'LLKA_tracepointInfo': (LLKA_TracepointInfos, []),
    'LLKA_tracepointState': (LLKA_Bool, [c_int32]),
}

_LIBRARY_NAMES = {
    'win32': ('LLKA.dll', 'libLLKA.dll'),
    'darwin': ('libLLKA.dylib',),
}


def _candidate_paths() -> list[Path]:
    # An explicitly configured library must not be silently replaced by another one
    env = os.environ.get('LLKA_LIBRARY')
    if env:
        return [Path(env)]

    candidates: list[Path] = []
    pkg_dir = Path(__file__).resolve().parent
    for name in _LIBRARY_NAMES.get(sys.platform, ('libLLKA.so',)):
        candidates.append(pkg_dir / name)
    # Versioned shared objects, e.g. libLLKA.so.0.1
    candidates.extend(sorted(pkg_dir.glob('libLLKA.so.*')))

    found = ctypes.util.find_library('LLKA')
    if found:
        candidates.append(Path(found))

    return candidates


def _load_library() -> ctypes.CDLL:
    # LLKA_CC is __stdcall on Windows which only matters on 32-bit x86
    loader = ctypes.WinDLL if sys.platform == 'win32' and struct.calcsize('P') == 4 else ctypes.CDLL  # type: ignore[attr-defined]

    tried: list[str] = []
    for path in _candidate_paths():
        if not (path.is_file() or path.name == str(path)):
            tried.append(str(path))
            continue
        try:
            return loader(str(path))
        except OSError as ex:
            tried.append(f'{path} ({ex})')

    raise ImportError(
        'Cannot load the libLLKA shared library. Build libLLKA with -DBUILD_SHARED_LIBRARY=ON and either '
        'install the Python package with the library or point the LLKA_LIBRARY environment variable to it. '
        'Tried: ' + ', '.join(tried)
    )


lib = _load_library()

for _name, (_restype, _argtypes) in PROTOTYPES.items():
    try:
        _func = getattr(lib, _name)
    except AttributeError:
        if _name in FILESYSTEM_FUNCTIONS:
            continue
        raise ImportError(f'The loaded libLLKA library does not export {_name}') from None
    _func.restype = _restype
    _func.argtypes = _argtypes


def has_function(name: str) -> bool:
    """Returns whether the loaded library exports a given C function."""
    return hasattr(lib, name)


def require_function(name: str):
    """Returns a C function or raises if this build of libLLKA does not provide it."""
    try:
        return getattr(lib, name)
    except AttributeError:
        from .errors import LLKANotImplementedError
        raise LLKANotImplementedError(
            f'{name} is not available in this build of libLLKA (filesystem access is disabled)'
        ) from None


def error_to_string(code: int) -> str:
    return lib.LLKA_errorToString(code).decode('ascii')


def check(code: int, context: str | None = None) -> None:
    """Raises an exception if ``code`` is not ``LLKA_OK``."""
    if code == 0:
        return

    message = error_to_string(code)
    if context:
        message = f'{context}: {message}'
    raise error_for_retcode(code, message)


def take_string(ptr: c_void_p | int | None) -> str | None:
    """Converts a string allocated by libLLKA to Python string and frees the original."""
    address = ptr.value if isinstance(ptr, c_void_p) else ptr
    if not address:
        return None
    try:
        return ctypes.string_at(address).decode('utf-8', errors='replace')
    finally:
        lib.LLKA_destroyString(address)


def encode(s: str) -> bytes:
    return s.encode('utf-8')


def decode(b: bytes | None) -> str:
    return '' if b is None else b.decode('utf-8', errors='replace')


def path_arg(path: str | os.PathLike[str]) -> bytes | str:
    """Converts a path to the representation expected by LLKA_PathChar parameters."""
    if sys.platform == 'win32':
        return os.fspath(path)
    return os.fsencode(path)

