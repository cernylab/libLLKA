Python examples
===

Python counterparts of the console examples in [examples](../../examples), written
with the `llka` Python bindings. Each program produces the same output as the
corresponding C program.

- [mirror_cif.py](mirror_cif.py)

  Loads and parses a mmCIF file and prints its contents back to standard output.

- [simple_ntc_assignment.py](simple_ntc_assignment.py)

  Reads a structure from a mmCIF file, breaks the structure down to __dinucleotide steps__,
  performs NtC assignment for each step and prints the results as a simple list.

- [classify_and_write_cif.py](classify_and_write_cif.py)

  Performs the same operation as `simple_ntc_assignment.py` but prints the input mmCIF data
  extended with categories that contain the results of the NtC assignment.

[dnatco_step_name.py](dnatco_step_name.py) is a helper module that constructs DNATCO step
names, the counterpart of `dnatco_step_name.c` used by the C examples.

Running the examples
---
The `llka` package must be importable, see [python/README.md](../README.md). For example,
from a build directory configured with `-DBUILD_PYTHON_BINDINGS=ON`:

    export PYTHONPATH=/path/to/libLLKA/python LLKA_LIBRARY=/path/to/build/libLLKA.so
    python3 /path/to/libLLKA/python/examples/simple_ntc_assignment.py structure.cif

The NtC assignment examples need the [assignment parametrization files](../../README.md#NtC-assignment-parametrization)
in the working directory. `classify_and_write_cif.py` also looks for them in `${CCP4}/share/dnatco`
and in the directory set by the `DNATCO_ASSETS_PATH` environment variable.
