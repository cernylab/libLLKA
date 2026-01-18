cmake_minimum_required(VERSION 3.15)

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG QUIET)

if (NOT pybind11_FOUND)
    # Try to find pybind11 via pip install location
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE pybind11_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if (pybind11_DIR)
        find_package(pybind11 CONFIG REQUIRED PATHS "${pybind11_DIR}")
    else ()
        message(FATAL_ERROR "pybind11 not found. Install with: pip install pybind11")
    endif ()
endif ()

# Python module target
pybind11_add_module(_pyllka_core
    src/bindings.cpp
)

# Link against libLLKA static library
target_link_libraries(_pyllka_core PRIVATE libLLKA_STATIC)

# Add compiler definitions (same as main library)
target_compile_definitions(_pyllka_core PRIVATE -DLLKA_COMPILER_GCC_LIKE)

# Include directories
target_include_directories(_pyllka_core PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}
    ${EIGEN_INCLUDE_DIR}
)

# Set output directory
set_target_properties(_pyllka_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pyllka"
)

# Copy Python package files to build directory
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/pyllka/__init__.py"
    "${CMAKE_CURRENT_BINARY_DIR}/pyllka/__init__.py"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/pyllka/_version.py"
    "${CMAKE_CURRENT_BINARY_DIR}/pyllka/_version.py"
    COPYONLY
)

# Copy asset files
file(COPY "${CMAKE_SOURCE_DIR}/assets/"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pyllka/assets/"
     FILES_MATCHING PATTERN "*.csv")

message(STATUS "Python bindings will be built: ${CMAKE_CURRENT_BINARY_DIR}/pyllka")
message(STATUS "Python executable: ${Python_EXECUTABLE}")
message(STATUS "pybind11 version: ${pybind11_VERSION}")
