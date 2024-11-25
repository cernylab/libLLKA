cmake_minimum_required(VERSION 3.20)

project(libLLKA LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(EIGEN_INCLUDE_DIR "" CACHE PATH "Path to the Eigen library files")
set(BUILD_STATIC_LIBRARY ON CACHE BOOL "Build static library")
set(BUILD_SHARED_LIBRARY OFF CACHE BOOL "Build shared library")
set(USE_ADDRESS_SANITIZER OFF CACHE BOOL "Use AddressSanitizer (supported only with recent Clang and GCC)")
set(BUILD_TESTING ON CACHE BOOL "Build tests")
set(BUILD_EXAMPLES ON CACHE BOOL "Build examples")
if (EMSCRIPTEN)
    set(
        EMX_JS_BUILD_MODE
        "ES6"
        CACHE STRING
        "How to build thd JavaScript module. Use 'Node' or 'Node_SIMD' for module that can be directly run in node, 'ES6' to build an ES6 module, 'ES6_NOSIMD' to build the ES6 module without SIMD code or 'ES6_NOWASM' to build the ES6 module as plain JavaScript code"
    )
    set(
        DNATCO_EXTRAS
        ON
        CACHE BOOL
        "Additional helper functions specific for the DNATCO tool"
    )
    set(
        ENABLE_WASM_PROFILING
        OFF
        CACHE BOOL
        "Generate DRAWF debugging information to enable profiling"
    )
else ()
    set(
        DNATCO_EXTRAS
        OFF
        CACHE BOOL
        "Additional helper functions specific for the DNATCO tool"
    )
endif ()

set(LIBLLKA_GLOBAL_DEFINITIONS -D_USE_MATH_DEFINES)
set(LIBLLKA_BUILD_DEFINITIONS -DLLKA_DLL_BUILD)
set(LIBLLKA_PLATFORM_DEFINITIONS "")

if (EMSCRIPTEN)
    set(BUILD_STATIC_LIBRARY ON)
    set(BUILD_SHARED_LIBRARY OFF)

    set("LLKA_PLATFORM_EMSCRIPTEN" "1")
    set(LIBLLKA_PLATFORM_DEFINITIONS -DLLKA_COMPILER_GCC_LIKE)

    if (EMX_JS_BUILD_MODE STREQUAL "Node")
        set(EMX_EMCC_WASM 1)
        set(EMX_EMCC_MODULARIZE 1)
        set(EMX_EMCC_ES6 0)
        set(EMX_EMCC_SIMD 0)
        set(EMX_EMCC_ENV node)
        set(EMX_EMCC_SUFFIX js)
        set(EMX_EMCC_NODERAWFS 1)
        add_link_options(-lnodefs.js -lnoderawfs.js)
    elseif (EMX_JS_BUILD_MODE STREQUAL "Node_SIMD")
        set(EMX_EMCC_WASM 1)
        set(EMX_EMCC_MODULARIZE 1)
        set(EMX_EMCC_ES6 0)
        set(EMX_EMCC_SIMD 1)
        set(EMX_EMCC_ENV node)
        set(EMX_EMCC_SUFFIX js)
        set(EMX_EMCC_NODERAWFS 1)
        add_definitions(-DLLKA_USE_SIMD_WASM)
        add_link_options(-lnodefs.js -lnoderawfs.js)
    elseif (EMX_JS_BUILD_MODE STREQUAL "ES6")
        set(EMX_EMCC_WASM 1)
        set(EMX_EMCC_MODULARIZE 1)
        set(EMX_EMCC_ES6 1)
        set(EMX_EMCC_SIMD 1)
        set(EMX_EMCC_ENV web)
        set(EMX_EMCC_SUFFIX js)
        set(EMX_EMCC_NODERAWFS 0)
        add_definitions(-DLLKA_FILESYSTEM_ACCESS_DISABLED)
    elseif (EMX_JS_BUILD_MODE STREQUAL "ES6_NOWASM")
        set(EMX_EMCC_WASM 0)
        set(EMX_EMCC_MODULARIZE 1)
        set(EMX_EMCC_ES6 1)
        set(EMX_EMCC_SIMD 0)
        set(EMX_EMCC_ENV web)
        set(EMX_EMCC_SUFFIX js)
        set(EMX_EMCC_NODERAWFS 0)
        add_definitions(-DLLKA_FILESYSTEM_ACCESS_DISABLED)
    elseif (EMX_JS_BUILD_MODE STREQUAL "ES6_NOSIMD")
        set(EMX_EMCC_WASM 1)
        set(EMX_EMCC_MODULARIZE 1)
        set(EMX_EMCC_ES6 1)
        set(EMX_EMCC_SIMD 0)
        set(EMX_EMCC_ENV web)
        set(EMX_EMCC_SUFFIX js)
        set(EMX_EMCC_NODERAWFS 0)
        add_definitions(-DLLKA_USE_SIMD_WASM)
        add_definitions(-DLLKA_FILESYSTEM_ACCESS_DISABLED)
    else()
        message(FATAL_ERROR "EMX_JS_BUILD_MODE ${EMX_JS_BUILD_MODE} is an invalid Emscripten build mode")
    endif ()

    # I honestly do not know if NO_DISABLE_EXCEPTION_CATCHING must passed both to the compiler and binding generator.
    # Passing it to both makes C++ try {} catch {} work. You may want to revise this at some point.
    add_compile_options(-fexceptions -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -fcolor-diagnostics)

    if (EMX_EMCC_SIMD EQUAL 1)
        add_compile_options(-msimd128)
    endif ()
elseif (WIN32)
    set("LLKA_PLATFORM_WIN32" "1")

    set(LIBLLKA_PLATFORM_DEFINITIONS "${LIBLLKA_PLATFORM_DEFINITIONS}" -DWIN32_LEAN_AND_MEAN -DNOMINMAX)

    if (MINGW OR MSYS)
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0.0)
            message(FATAL_ERROR "MinGW 7.0.0 or newer is required to build libLLKA")
        endif ()
        add_compile_options(-Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -g -msse2)
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            add_compile_options(-O2)
        endif ()

        if (USE_ADDRESS_SANITIZER)
            add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize=leak -fstack-protector-strong)
            add_link_options(-fsanitize=address -fsanitize=undefined -fsanitize=leak)
        endif ()

        if (MINGW)
            set(LIBLLKA_PLATFORM_DEFINITIONS -DLLKA_COMPILER_MINGW)
        else ()
            set(LIBLLKA_PLATFORM_DEFINITIONS -DLLKA_COMPILER_MSYS)
        endif ()

        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (MSVC)
        set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;MinSizeRel")

        set(CMAKE_CXX_FLAGS_RELEASE "/GL")
        set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "/INCREMENTAL:NO /LTCG")
        set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "/INCREMENTAL:NO /LTCG")

        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/GL")
        set(CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO "/debug /INCREMENTAL:NO /LTCG")
        set(CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO "/debug /INCREMENTAL:NO /LTCG")
        set(LIBLLKA_PLATFORM_DEFINITIONS "${LIBLLKA_PLATFORM_DEFINITIONS}" -DLLKA_COMPILER_MSVC -DLLKA_HAVE_SAD_LADY_COMPILER -DUNICODE -D_UNICODE)
        set(DEFAULT_SYMVER_LINK " ")
    else()
        message(FATAL_ERROR "Compiler not recognized or not supported")
    endif()
elseif (UNIX)
    set("LLKA_PLATFORM_UNIX" "1")

    if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        add_compile_options(-msse2)
        add_definitions("-DLLKA_USE_SIMD_X86")
    endif ()

    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        add_compile_options(-pipe -Wall -Wextra -pedantic -fvisibility=hidden -g)
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND
                CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.1.0 OR
                CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 9.1.0)
                add_compile_options(-flto)
            endif ()
            add_compile_options(-O2)
        elseif (CMAKE_BUILD_TYPE STREQUAL "Profile")
            add_compile_options(-O2 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer)
        endif ()

        if (USE_ADDRESS_SANITIZER)
            add_compile_options(-fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize=leak -fstack-protector-strong)
            add_link_options(-fsanitize=address -fsanitize=undefined -fsanitize=leak)
        endif ()

        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
            add_compile_options(-fdiagnostics-color=always)
        elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
            add_compile_options(-fcolor-diagnostics)
        endif ()

        set(LIBLLKA_PLATFORM_DEFINITIONS -DLLKA_COMPILER_GCC_LIKE)
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    endif ()
else()
    message(FATAL_ERROR "Platform not recognized or not supported")
endif()

if ((NOT BUILD_STATIC_LIBRARY) AND (NOT BUILD_SHARED_LIBRARY))
    message(FATAL "Neither static nor shared library is set to be built. There is nothing to do.")
endif ()

include_directories(
    ${INCLUDE_DIRECTORIES}
    include
    src
    "${CMAKE_CURRENT_BINARY_DIR}"
    SYSTEM ${EIGEN_INCLUDE_DIR}
)

set(
    libLLKA_SRCS
    "src/minicif/parser.cpp"
    "src/minicif/writer.cpp"
    "src/util/elementaries.cpp"
    "src/util/geometry.cpp"
    "src/util/geometry_cpp.cpp"
    "src/util/printers.cpp"
    "src/classification.cpp"
    "src/connectivity_similarity.cpp"
    "src/extend.cpp"
    "src/extract.cpp"
    "src/measurements.cpp"
    "src/minicif.cpp"
    "src/ntc.cpp"
    "src/ntc_bones.cpp"
    "src/ntc_constants.cpp"
    "src/nucleotide.cpp"
    "src/llka.cpp"
    "src/resource_loaders.cpp"
    "src/segmentation.cpp"
    "src/structure.cpp"
    "src/structure_util.cpp"
    "src/superposition.cpp"
    "src/tracing.cpp"
    "src/llka_cpp.cpp"
)

set(
    libLLKA_PUB_HDRS
    "include/llka_classification.h"
    "include/llka_cpp.h"
    "include/llka_connectivity_similarity.h"
    "include/llka_main.h"
    "include/llka_minicif.h"
    "include/llka_measurements.h"
    "include/llka_module.h"
    "include/llka_ntc.h"
    "include/llka_nucleotide.h"
    "include/llka_resource_loaders.h"
    "include/llka_segmentation.h"
    "include/llka_structure.h"
    "include/llka_superposition.h"
    "include/llka_tracing.h"
    "include/llka_util.h"
)

if (DNATCO_EXTRAS)
    set(
        libLLKA_SRCS
        ${libLLKA_SRCS}
        "dnatco_extras/src/dnatco.cpp"
    )
endif ()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/llka_config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/llka_config.h")

set(LLKA_TARGETS "")

if (EMSCRIPTEN)
    if (EMX_JS_BUILD_MODE STREQUAL "ES6_NOWASM")
        set(WASM_SUFFIX "_nowasm")
        set(TARGET_OUTPUT_NAME LLKA_nowasm)
    elseif (EMX_JS_BUILD_MODE STREQUAL "ES6_NOSIMD")
        set(WASM_SUFFIX "_nosimd")
        set(TARGET_OUTPUT_NAME LLKA_nosimd)
    elseif ((EMX_JS_BUILD_MODE STREQUAL "Node") OR (EMX_JS_BUILD_MODE STREQUAL "Node_SIMD"))
        set(WASM_SUFFIX "_node")
        set(TARGET_OUTPUT_NAME LLKA_node)
    else ()
        set(WASM_SUFFIX "")
        set(TARGET_OUTPUT_NAME LLKA)
    endif ()

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
        set(EMX_OPTIFLAGS "-g -Og")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -flto")
        set(EMX_OPTIFLAGS "-O3")
        set(EMX_OPTIFLAGS_LTO "-flto")
    endif ()

    if (ENABLE_WASM_PROFILING)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
        set(EMX_OPTIFLAGS "-g")
    endif ()

    add_library(libLLKA_STATIC STATIC ${libLLKA_SRCS})
    set_target_properties(
        libLLKA_STATIC
        PROPERTIES
            OUTPUT_NAME ${TARGET_OUTPUT_NAME}
            LINK_OPTIONS "--bind -fexceptions -sNO_DISABLE_EXCEPTION_CATCHING -sNODERAWFS=${EMX_EMCC_NODERAWFS} -sENVIRONMENT=${EMX_EMCC_ENV} -sWASM=${EMX_EMCC_WASM} -sALLOW_MEMORY_GROWTH -sMODULARIZE=${EMX_EMCC_MODULARIZE} -sEXPORT_ES6=${EMX_EMCC_ES6} ${EMX_OPTIFLAGS} ${EMX_OPTIFLAGS_LTO}"
    )
    target_compile_definitions(libLLKA_STATIC PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_BUILD_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
    target_link_libraries(libLLKA_STATIC ${LLKA_EXTRA_LINK_LIBS})

    set(LLKA_TARGETS ${LLKA_TARGETS} libLLKA_STATIC)

    add_custom_command(
        TARGET libLLKA_STATIC POST_BUILD
        COMMAND em++ --bind -fexceptions -sNODERAWFS=${EMX_EMCC_NODERAWFS} -sENVIRONMENT=${EMX_EMCC_ENV} -sWASM=${EMX_EMCC_WASM} -sALLOW_MEMORY_GROWTH -sMODULARIZE=${EMX_EMCC_MODULARIZE} -sEXPORT_ES6=${EMX_EMCC_ES6} ${EMX_OPTIFLAGS} ${EMX_OPTIFLAGS_LTO} -o "${CMAKE_CURRENT_BINARY_DIR}/libLLKA${WASM_SUFFIX}.${EMX_EMCC_SUFFIX}" -Wl,--whole-archive "${CMAKE_CURRENT_BINARY_DIR}/libLLKA${WASM_SUFFIX}.a" -Wl,--no-whole-archive
        COMMENT "Generating JavaScript bindings"
    )
else ()
    if (BUILD_STATIC_LIBRARY)
        add_library(libLLKA_STATIC STATIC ${libLLKA_SRCS})
        set_target_properties(
            libLLKA_STATIC
            PROPERTIES OUTPUT_NAME LLKA
            LINK_FLAGS ${DEFAULT_SYMVER_LINK}
        )
        target_compile_definitions(libLLKA_STATIC PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_BUILD_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
        target_link_libraries(libLLKA_STATIC ${LLKA_EXTRA_LINK_LIBS})

        set(LLKA_TARGETS ${LLKA_TARGETS} libLLKA_STATIC)
    endif ()

    if (BUILD_SHARED_LIBRARY)
        add_library(libLLKA_SHARED SHARED ${libLLKA_SRCS})
        set_target_properties(
            libLLKA_SHARED
            PROPERTIES VERSION 0.1
                       SOVERSION 0.1
                       OUTPUT_NAME LLKA
                       LINK_FLAGS ${DEFAULT_SYMVER_LINK}
        )
        target_compile_definitions(libLLKA_SHARED PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_BUILD_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
        target_link_libraries(libLLKA_SHARED ${LLKA_EXTRA_LINK_LIBS})

        set(LLKA_TARGETS ${LLKA_TARGETS} libLLKA_SHARED)
    endif ()
endif ()

if (UNIX)
    include(GNUInstallDirs)
else ()
    set(CMAKE_INSTALL_LIBDIR "lib")
    set(CMAKE_INSTALL_INCLUDEDIR "include")
    set(CMAKE_INSTALL_BINDIR "bin")
endif ()

install(
    TARGETS ${LLKA_TARGETS} EXPORT libLLKA_targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(
    EXPORT libLLKA_targets DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/LLKA"
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/llka_config.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/LLKA"
)
install(
    FILES ${libLLKA_PUB_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/LLKA"
)
if (DNATCO_EXTRAS)
    install(
        FILES "${CMAKE_CURRENT_SOURCE_DIR}/dnatco_extras/include/dnatco.h"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/LLKA/"
    )
endif ()

if (EMSCRIPTEN)
    install(
        FILES "${CMAKE_CURRENT_BINARY_DIR}/libLLKA${WASM_SUFFIX}.${EMX_EMCC_SUFFIX}"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    )
    if (EMX_EMCC_WASM EQUAL 1)
        install(
            FILES "${CMAKE_CURRENT_BINARY_DIR}/libLLKA${WASM_SUFFIX}.wasm"
            DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        )
    endif ()
endif ()

if (BUILD_TESTING)
    include(CTest)
    add_subdirectory("tests")
endif ()

if (BUILD_EXAMPLES)
    add_subdirectory("examples")
endif ()

if (WIN32)
    set(LIBLLKA_CMAKE_DIR cmake)
else ()
    set(LIBLLKA_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/LLKA/")
endif ()

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/LLKAConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/LLKAConfigForInstall.cmake"
    INSTALL_DESTINATION "${LIBLLKA_CMAKE_DIR}"
    PATH_VARS
        CMAKE_INSTALL_PREFIX
        CMAKE_INSTALL_INCLUDEDIR
        CMAKE_INSTALL_LIBDIR
        CMAKE_INSTALL_BINDIR
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/LLKAConfigForInstall.cmake"
    DESTINATION "${LIBLLKA_CMAKE_DIR}"
    RENAME LLKAConfig.cmake
)

if (WIN32)
elseif (APPLE)
else ()
    set(PKGCONFIG_PLATFORM_LIBS "-lm")
endif ()

configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/pkgconfig/llka.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/llka.pc" @ONLY
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/llka.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/"
)
