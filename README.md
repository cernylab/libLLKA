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
