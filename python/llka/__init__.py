"""Python bindings for libLLKA, the (Re)DNATCO nucleic acid structure processing library.

The bindings call the libLLKA C API through :mod:`ctypes`. Errors reported by
libLLKA are raised as exceptions derived from :class:`LLKAError`.

The submodules mirror the C headers:

- :mod:`llka.core` - points and matrices (``llka_main.h``)
- :mod:`llka.structure` - atoms and structures (``llka_structure.h``)
- :mod:`llka.measurements` - angles and distances between atoms (``llka_measurements.h``)
- :mod:`llka.ntc` - NtC steps and their metrics (``llka_ntc.h``)
- :mod:`llka.nucleotide` - nucleotides and sugar puckers (``llka_nucleotide.h``)
- :mod:`llka.superposition` - superposition and RMSD (``llka_superposition.h``)
- :mod:`llka.segmentation` - models, chains and residues (``llka_segmentation.h``)
- :mod:`llka.connectivity_similarity` - step connectivity and similarity (``llka_connectivity_similarity.h``)
- :mod:`llka.classification` - NtC classification (``llka_classification.h``)
- :mod:`llka.resources` - loading of classification parameters (``llka_resource_loaders.h``)
- :mod:`llka.minicif` - mmCIF input and output (``llka_minicif.h``)
- :mod:`llka.tracing` - debug tracing (``llka_tracing.h``)

All public names are also available directly from the ``llka`` package.
"""

from __future__ import annotations

from . import (
    classification,
    connectivity_similarity,
    core,
    enums,
    errors,
    measurements,
    minicif,
    ntc,
    nucleotide,
    resources,
    segmentation,
    structure,
    superposition,
    tracing,
)
from .classification import *  # noqa: F403
from .connectivity_similarity import *  # noqa: F403
from .core import *  # noqa: F403
from .enums import *  # noqa: F403
from .errors import *  # noqa: F403
from .measurements import *  # noqa: F403
from .minicif import *  # noqa: F403
from .ntc import *  # noqa: F403
from .nucleotide import *  # noqa: F403
from .resources import *  # noqa: F403
from .segmentation import *  # noqa: F403
from .structure import *  # noqa: F403
from .superposition import *  # noqa: F403
from .tracing import *  # noqa: F403

__version__ = '0.1.0'

__all__ = [
    '__version__',
    *classification.__all__,
    *connectivity_similarity.__all__,
    *core.__all__,
    *enums.__all__,
    *errors.__all__,
    *measurements.__all__,
    *minicif.__all__,
    *ntc.__all__,
    *nucleotide.__all__,
    *resources.__all__,
    *segmentation.__all__,
    *structure.__all__,
    *superposition.__all__,
    *tracing.__all__,
]
