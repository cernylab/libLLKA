#
# Python bindings for libLLKA
#
# The bindings are a pure Python package that loads the shared libLLKA library
# through ctypes. Nothing is compiled here; the package and the library are
# installed together and the Python test suite is registered with CTest.
#

if (NOT BUILD_SHARED_LIBRARY)
    message(
        FATAL_ERROR
        "Python bindings need the shared libLLKA library. Configure with -DBUILD_STATIC_LIBRARY=OFF -DBUILD_SHARED_LIBRARY=ON"
    )
endif ()

find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)

set(
    LLKA_PYTHON_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages"
    CACHE PATH
    "Directory where to install the llka Python package"
)

install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/llka"
    DESTINATION "${LLKA_PYTHON_INSTALL_DIR}"
    FILES_MATCHING
        PATTERN "*.py"
        PATTERN "py.typed"
        PATTERN "__pycache__" EXCLUDE
)
# The package looks for the library in its own directory first
install(
    FILES "$<TARGET_FILE:libLLKA_SHARED>"
    DESTINATION "${LLKA_PYTHON_INSTALL_DIR}/llka"
)

if (BUILD_TESTING)
    set(
        LLKA_PYTHON_TESTS
        test_abi
        test_api_coverage
        test_classification
        test_connectivity
        test_examples
        test_extract_backbone
        test_measurements
        test_minicif
        test_ntc
        test_nucleotide
        test_segmentation
        test_similarity
        test_structure_manipulation
        test_structure_splitting
        test_superposition
    )

    set(LLKA_PYTHON_TEST_ENV
        "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}"
        "LLKA_LIBRARY=$<TARGET_FILE:libLLKA_SHARED>"
        "PYTHONDONTWRITEBYTECODE=1"
    )

    # The ABI test compiles a small C program that checks the ctypes structure layouts
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
        string(JOIN " " LLKA_PYTHON_TEST_C_FLAGS ${LIBLLKA_PLATFORM_DEFINITIONS} "-I${PROJECT_BINARY_DIR}")
        list(APPEND LLKA_PYTHON_TEST_ENV
            "LLKA_TEST_C_COMPILER=${CMAKE_C_COMPILER}"
            "LLKA_TEST_C_FLAGS=${LLKA_PYTHON_TEST_C_FLAGS}"
        )
    endif ()

    # test_examples compares the output of the Python examples with the output of the C examples
    if (BUILD_EXAMPLES AND NOT EMSCRIPTEN)
        list(APPEND LLKA_PYTHON_TEST_ENV "LLKA_C_EXAMPLES_DIR=${PROJECT_BINARY_DIR}/examples")
    endif ()

    foreach (test_module ${LLKA_PYTHON_TESTS})
        add_test(
            NAME "Python_${test_module}"
            COMMAND "${Python3_EXECUTABLE}" -m unittest -v ${test_module}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        )
        set_tests_properties(
            "Python_${test_module}"
            PROPERTIES ENVIRONMENT "${LLKA_PYTHON_TEST_ENV}"
        )
    endforeach ()
endif ()
