"""
pyllka - Python bindings for libLLKA

A nucleic acid structure analysis library from the (Re)DNATCO project.
"""

from ._version import __version__
from ._pyllka_core import (
    # Core types
    Point,
    Atom,
    Structure,
    ImportedStructure,

    # Enums
    NtC,
    CANA,
    SugarPucker,
    SugarPuckerNameBrevity,
    RetCode,

    # Classification
    ClassificationContext,
    ClassifiedStep,
    StepMetrics,
    ConfalScore,
    AverageConfal,
    NuAngles,
    average_confal_attempted,

    # Similarity and Connectivity
    Similarity,
    Connectivity,
    measure_step_similarity_ntc,
    measure_step_similarity_ntc_multiple,
    measure_step_similarity_structure,
    measure_step_connectivity_ntcs,
    measure_step_connectivity_ntcs_multiple_first,
    measure_step_connectivity_ntcs_multiple_second,
    measure_step_connectivity_structures,

    # Exceptions
    LLKAError,

    # Utility functions
    ntc_to_name,
    name_to_ntc,
    cana_to_name,
    sugar_pucker_to_name,
    classification_violation_to_name,
    rad2deg,
    full_angle_from_rad,

    # CIF I/O
    load_structure_with_cif,
    cif_to_string,

    # Version constants
    NTC_VERSION,
    CANA_VERSION,
)

from .io import (
    load_structure,
    save_structure,
)

from .classification import (
    create_classification_context,
    classify_structure,
    classify_step,
)

from .utils import (
    update_coordinates,
    split_to_dinucleotides,
)

# Optional MDAnalysis integration
try:
    from .mdanalysis_integration import (
        TrajectoryAnalyzer,
        analyze_trajectory,
    )
    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False

__all__ = [
    '__version__',
    # Core types
    'Point',
    'Atom',
    'Structure',
    'ImportedStructure',
    # Enums
    'NtC',
    'CANA',
    'SugarPucker',
    'SugarPuckerNameBrevity',
    'RetCode',
    # Classification
    'ClassificationContext',
    'ClassifiedStep',
    'StepMetrics',
    'ConfalScore',
    'AverageConfal',
    'NuAngles',
    'create_classification_context',
    'classify_structure',
    'classify_step',
    'average_confal_attempted',
    # Similarity and Connectivity
    'Similarity',
    'Connectivity',
    'measure_step_similarity_ntc',
    'measure_step_similarity_ntc_multiple',
    'measure_step_similarity_structure',
    'measure_step_connectivity_ntcs',
    'measure_step_connectivity_ntcs_multiple_first',
    'measure_step_connectivity_ntcs_multiple_second',
    'measure_step_connectivity_structures',
    # I/O
    'load_structure',
    'save_structure',
    'load_structure_with_cif',
    'cif_to_string',
    # Utils
    'update_coordinates',
    'split_to_dinucleotides',
    'ntc_to_name',
    'name_to_ntc',
    'cana_to_name',
    'sugar_pucker_to_name',
    'classification_violation_to_name',
    'rad2deg',
    'full_angle_from_rad',
    # Version constants
    'NTC_VERSION',
    'CANA_VERSION',
    # Classification violation constants
    'CLASSIFICATION_OK',
    'CLASSIFICATION_E_SCORE_TOO_LOW',
    'CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS',
    'CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT',
    'CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT',
    'CLASSIFICATION_E_CC_TOO_LOW',
    'CLASSIFICATION_E_CC_TOO_HIGH',
    'CLASSIFICATION_E_NN_TOO_LOW',
    'CLASSIFICATION_E_NN_TOO_HIGH',
    'CLASSIFICATION_E_MU_TOO_LOW',
    'CLASSIFICATION_E_MU_TOO_HIGH',
    'CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH',
    'CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT',
    'CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT',
    'CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES',
    'CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED',
    'CLASSIFICATION_E_WRONG_METRICS',
    'CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH',
    # Exceptions
    'LLKAError',
]

if _HAS_MDANALYSIS:
    __all__.extend([
        'TrajectoryAnalyzer',
        'analyze_trajectory',
    ])

# Import classification violation constants (at end to avoid circular import)
# Access the already-imported _pyllka_core module
import sys
_core = sys.modules.get('pyllka._pyllka_core')
if _core is not None:
    CLASSIFICATION_OK = getattr(_core, 'CLASSIFICATION_OK', 0)
    CLASSIFICATION_E_SCORE_TOO_LOW = getattr(_core, 'CLASSIFICATION_E_SCORE_TOO_LOW', 1)
    CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS = getattr(_core, 'CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS', 2)
    CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT = getattr(_core, 'CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT', 4)
    CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT = getattr(_core, 'CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT', 8)
    CLASSIFICATION_E_CC_TOO_LOW = getattr(_core, 'CLASSIFICATION_E_CC_TOO_LOW', 16)
    CLASSIFICATION_E_CC_TOO_HIGH = getattr(_core, 'CLASSIFICATION_E_CC_TOO_HIGH', 32)
    CLASSIFICATION_E_NN_TOO_LOW = getattr(_core, 'CLASSIFICATION_E_NN_TOO_LOW', 64)
    CLASSIFICATION_E_NN_TOO_HIGH = getattr(_core, 'CLASSIFICATION_E_NN_TOO_HIGH', 128)
    CLASSIFICATION_E_MU_TOO_LOW = getattr(_core, 'CLASSIFICATION_E_MU_TOO_LOW', 256)
    CLASSIFICATION_E_MU_TOO_HIGH = getattr(_core, 'CLASSIFICATION_E_MU_TOO_HIGH', 512)
    CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH = getattr(_core, 'CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH', 1024)
    CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT = getattr(_core, 'CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT', 2048)
    CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT = getattr(_core, 'CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT', 4096)
    CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES = getattr(_core, 'CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES', 8192)
    CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED = getattr(_core, 'CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED', 16384)
    CLASSIFICATION_E_WRONG_METRICS = getattr(_core, 'CLASSIFICATION_E_WRONG_METRICS', 32768)
    CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH = getattr(_core, 'CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH', 65536)
