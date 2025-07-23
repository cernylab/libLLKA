libLLKA - The (Re)DNATCO nucleic acid processing library
===

Dependencies
---

libLLKA depends on the following 3rd-party libraries

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library. libLLKA is known to work with Eigen version 3.4.0.

How to build
---
libLLKA uses [CMake](https://cmake.org/) build system to build itself.

On a UNIX system, issue the following commands to build libLLKA

    cd <path_to_libLLKA_source_code_directory>
    mkdir build
    cd build
    cmake .. -DEIGEN_INCLUDE_DIR=/path/to/Eigen
    make

Compile options
---
- `-DBUILD_STATIC_LIBRARY` `[ON|OFF]` Whether to build or do not build a binary that can be linked statically
- `-DBUILD_SHARED_LIBRARY` `[ON|OFF]` Whether to build or do not build a binary that can be linked dynamically
- `-DUSE_ADDRESS_SANITIZER` `[ON|OFF]` Whether to use [ASan](https://github.com/google/sanitizers/wiki/AddressSanitizer) to check for all kinds of memory access errors and undefined behavior. May not be available on all platforms.

- `-DBUILD_EXAMPLES` `[ON|OFF]` Whether to build examples. Note that some examples require additional dependencies to build. Consult the README files in the directories of the individual examples for more details.

Compiling with Emscripten
---
Emscripten is rather quirky when it comes to dealing with const vs. non-const member functions, classes with explicitly deleted constructors or assignment operators and other
"advanced" C++ techniques. libLLKA is known to build and run when compiled with Emscripten 3.1.59. Use of other versions of Emscripten may result in build failures.

NtC assignment parametrization
---
The NtC assignment process uses a series of parameters whose values affect the results. libLLKA needs to load these parameters before any NtC assignment can be performed. Currently, the parameters are defined in 5 CSV files that can be split to two categories:

#### Parameters used to determine the NtC and CANA class of a __dinucleotide step__
- [golden_steps.csv](assets/golden_steps.csv)
- [clusters.csv](assets/clusters.csv)
- [nu_angles.csv](assets/nu_angles.csv)

#### Parameters used to calculate Confal score for a __dinucleotide step__
- [confals.csv](assets/confals.csv)
- [confal_percentiles.csv](assets/confal_percentiles.csv)

Note that any program that uses libLLKA to perform the NtC assignment will need to load the assignment parameters and, in turn, will likely need to access these files.

Example code and tools
---
Examples provided in the [examples](examples/) directory are primarily intended for programmers who would like to use libLLKA in their programs. The examples provide a basic overview of how to load a structure, perform the NtC assignemt, retrieve the results and handle various errors that may occur throughout the process.

There are two example programs that may be used as a very simplifed version of the [DNATCO web-based tool](https://dnatco.datmos.org).

 - [classify_and_write_cif](examples/classify_and_write_cif)

   A CLI-only tool that reads a structure from a mmCIF file, performs the NtC assignment and writes the results into a mmCIF file with additional categories that contain the NtC assignment results.

 - [similarity_connectivity](examples/similarity_connectivity)

   A CLI-only tool that reads a structure from a mmCIF file, performs the NtC assignment and writes the results into a mmCIF file with additional categories that contain the NtC assignment results. It also produces the "similarity" and "connectivity" data in JSON format. The level of details can be controled by additional commandline switches, see --help.

 - [gui_assigner](examples/gui_assigner)

   A graphical tool that has a similar functionality but also provides a GUI to enter data and display results. Results are displayed in a table, double-clicking on a row reveals more details about the assigned step.

   The `gui_assigner` example also serves a secondary purpose of demonstrating the use of libLLKAs's C++ interface.

   `gui_assigner` relies on [Qt 6](https://www.qt.io/product/qt6) libraries. Qt 6 runtime and development files must be present on the system to build and run the `gui_assigner` example.

Intentionally unimplemented functionality
---
libLLKA is written to be a standalone, self-sufficient library that can be easily integrated into tools for structural analysis of biomolecules. It implements only the functionality that is essential for its purpose and limits dependencies on 3rd-party code as much as it is practically possible. The __notably unsupported__ features include:

  - Reading of structures from PDB files. The PDB file format is now widely considered a legacy format that should not be used for new structures. It is assumed that a tool using libLLKA will be able to convert from PDB to mmCIF or create libLLKA structure object directly from its own internal representation.
  - Reading of compressed files.
