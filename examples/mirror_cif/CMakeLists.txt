if (BUILD_STATIC_LIBRARY)
    set(LLKA_LIB_LINK libLLKA_STATIC)
else ()
    set(LLKA_LIB_LINK libLLKA_SHARED)
endif ()

# Simple NtC assignment example
add_executable(
    mirror_cif
    mirror_cif.c
)
set_target_properties(
    mirror_cif
    PROPERTIES
        LANGUAGE C
)
target_compile_definitions(mirror_cif PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})

if (EMSCRIPTEN)
    set(MAYBE_EMX_LINK_FLAGS "-sALLOW_MEMORY_GROWTH -fexceptions")
    set(MAYBE_EMX_LINK_LIBS ${LLKA_LIB_LINK} nodefs.js noderawfs.js)
else ()
    set(MAYBE_EMX_LINK_FLAGS "")
    set(MAYBE_EMX_LINK_LIBS "${LLKA_LIB_LINK}")
endif ()

set_target_properties(
    mirror_cif
    PROPERTIES
        LANGUAGE C
        LINK_FLAGS "${MAYBE_EMX_LINK_FLAGS}"
)
target_link_libraries(mirror_cif ${MAYBE_EMX_LINK_LIBS})
