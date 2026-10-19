llka - Python bindings for libLLKA
===

`llka` is a pure Python package that calls the libLLKA C API through
[ctypes](https://docs.python.org/3/library/ctypes.html). It needs Python 3.10
or newer and has no dependencies besides the Python standard library and the
libLLKA shared library.

Building and installing
---
The bindings need the shared variant of libLLKA:

    mkdir build
    cd build
    cmake .. -DEIGEN_INCLUDE_DIR=/path/to/Eigen \
             -DBUILD_STATIC_LIBRARY=OFF -DBUILD_SHARED_LIBRARY=ON \
             -DBUILD_PYTHON_BINDINGS=ON
    make
    make install

`make install` puts the `llka` package, together with a copy of the shared
library, into `<prefix>/<libdir>/pythonX.Y/site-packages`. Use
`-DLLKA_PYTHON_INSTALL_DIR=...` to install it elsewhere.

To use the package without installing it, add the `python` directory to
`PYTHONPATH` and point the `LLKA_LIBRARY` environment variable to the built
library:

    PYTHONPATH=/path/to/libLLKA/python LLKA_LIBRARY=/path/to/build/libLLKA.so python3

Running the tests
---
With `BUILD_TESTING` enabled, `ctest` runs the Python test suite along with the
C++ tests. The Python tests mirror the C++ tests in `tests/` and use the same
test data. `test_abi` additionally compiles a small C program to check that the
ctypes mirrors of the C structures and enums match the C headers, and
`test_api_coverage` checks that every callable function of the C API is used by
the bindings.

Example
---
Complete example programs are in [examples](examples). They are Python
counterparts of the console examples of libLLKA. A short example follows:

```python
import math
import llka

def load(name, kind):
    return llka.load_resource_file(f'assets/{name}', kind)

ctx = llka.ClassificationContext(
    clusters=load('clusters.csv', llka.ResourceType.CLUSTERS),
    golden_steps=load('golden_steps.csv', llka.ResourceType.GOLDEN_STEPS),
    confals=load('confals.csv', llka.ResourceType.CONFALS),
    cluster_nu_angles=load('nu_angles.csv', llka.ResourceType.AVERAGE_NU_ANGLES),
    confal_percentiles=load('confal_percentiles.csv', llka.ResourceType.CONFAL_PERCENTILES),
    limits=llka.ClassificationLimits(
        minimum_nearest_neighbors=7,
        number_of_used_nearest_neighbors=11,
        average_neighbors_torsion_cutoff=math.radians(28.0),
        nearest_neighbor_torsions_cutoff=math.radians(28.0),
        total_distance_cutoff=math.radians(60.0),
        pseudorotation_cutoff=math.radians(72.0),
        minimum_cluster_votes=0.001111,
    ),
    max_close_enough_rmsd=0.5,
)

imported = llka.cif_file_to_structure('assets/test_cifs/1BNA.cif')
steps = imported.structure.split_to_dinucleotide_steps()

for attempt in ctx.classify_steps(steps):
    if attempt.error is not None:
        print('Not a valid step:', attempt.error)
    elif attempt.step.is_assigned:
        print(attempt.step.assigned_ntc.name, attempt.step.confal_score.total)
    else:
        print('Unassigned, closest NtC', attempt.step.closest_ntc.name, attempt.step.violations)
```

Design
---
- **Modules mirror the C headers.** `llka.structure` wraps `llka_structure.h`,
  `llka.ntc` wraps `llka_ntc.h` and so on. All public names are also available
  directly from the `llka` package.
- **Errors are exceptions.** Every `LLKA_RetCode` other than `LLKA_OK` is raised
  as a subclass of `llka.LLKAError` that carries the code in its `retcode`
  attribute. The exceptions also derive from the fitting built-in exception,
  e.g. `ValueError` or `FileNotFoundError`. Sentinel values of the C API are
  translated as well: conversions of unknown names raise `ValueError` and a
  missing NtC or CANA in classification results is `None` instead of
  `NtC.INVALID` or `CANA.INVALID`.
- **Data stays in libLLKA memory.** `Structure`, `Points`, `Matrix`, `CifData`
  and `ResourceList` keep their data in memory managed by libLLKA and release it
  when they are garbage-collected. Atoms are converted to Python objects only
  when they are accessed, so large structures can be passed between libLLKA
  functions without conversion.
- **Small types are Python values.** `Atom`, `Point`, `StepMetrics`,
  `ClassifiedStep` and other plain C structs are immutable Python values.
  Changing an atom of a structure therefore means assigning a new atom:
  `stru[0] = stru[0].replace(label_seq_id=2)`.
- **Views are checked.** A `StructureView` refers to atoms of its parent
  `Structure` and keeps the parent alive. Adding or removing atoms in the parent
  invalidates its views; using an invalid view raises `RuntimeError` instead of
  touching freed memory. Cif categories and items are looked up by name on each
  access for the same reason.

C API functions that are not exposed
---
The `LLKA_MINICIF_NORMALIZE` import option is not exposed because it does not
order atoms correctly.

A few functions have no Python counterpart because the Python types make them
unnecessary, e.g. `LLKA_makeAtom` (an `Atom` is created directly). See
`UNWRAPPED_FUNCTIONS` in `llka/_capi.py`.
