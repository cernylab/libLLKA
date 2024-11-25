if (BUILD_STATIC_LIBRARY)
    set(LLKA_LIB_LINK libLLKA_STATIC)
else ()
    set(LLKA_LIB_LINK libLLKA_SHARED)
endif ()

# Simple NtC assignment example
add_executable(
    classify_and_write_cif
    classify_and_write_cif.c
    dnatco_step_name.c
)
target_compile_definitions(classify_and_write_cif PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})

if (EMSCRIPTEN)
    set(MAYBE_EMX_LINK_FLAGS "-sALLOW_MEMORY_GROWTH -fexceptions")
    set(MAYBE_EMX_LINK_LIBS ${LLKA_LIB_LINK} nodefs.js noderawfs.js)
else ()
    set(MAYBE_EMX_LINK_FLAGS "")
    set(MAYBE_EMX_LINK_LIBS "${LLKA_LIB_LINK}")
endif ()

set_target_properties(
    classify_and_write_cif
    PROPERTIES
        LANGUAGE C
        LINK_FLAGS "${MAYBE_EMX_LINK_FLAGS}"
)
target_link_libraries(classify_and_write_cif ${MAYBE_EMX_LINK_LIBS})

# Get definition files necessary to initialize classification context
add_custom_command(
    TARGET classify_and_write_cif POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/../../assets/clusters.csv"
    "${CMAKE_CURRENT_BINARY_DIR}/clusters.csv"
)
add_custom_command(
    TARGET classify_and_write_cif POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/../../assets/confals.csv"
    "${CMAKE_CURRENT_BINARY_DIR}/confals.csv"
)
add_custom_command(
    TARGET classify_and_write_cif POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/../../assets/golden_steps.csv"
    "${CMAKE_CURRENT_BINARY_DIR}/golden_steps.csv"
)
add_custom_command(
    TARGET classify_and_write_cif POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/../../assets/nu_angles.csv"
    "${CMAKE_CURRENT_BINARY_DIR}/nu_angles.csv"
)
add_custom_command(
    TARGET classify_and_write_cif POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/../../assets/confal_percentiles.csv"
    "${CMAKE_CURRENT_BINARY_DIR}/confal_percentiles.csv"
)
