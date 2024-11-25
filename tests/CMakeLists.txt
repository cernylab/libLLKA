if (BUILD_STATIC_LIBRARY)
    set(LLKA_LIB_LINK libLLKA_STATIC)
else ()
    set(LLKA_LIB_LINK libLLKA_SHARED)
endif ()

set(CPP_EXTRA_LINK "")
if (EMSCRIPTEN)
    set(CPP_EXTRA_LINK embind)
endif ()

add_executable(test_structure_manipulation test_structure_manipulation.cpp effedup.cpp)
target_compile_definitions(test_structure_manipulation PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_structure_manipulation ${LLKA_LIB_LINK})
add_test(NAME StructureManipulation COMMAND test_structure_manipulation)

add_executable(test_structure_splitting test_structure_splitting.cpp effedup.cpp)
target_compile_definitions(test_structure_splitting PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_structure_splitting PRIVATE ${LLKA_LIB_LINK})
add_test(NAME StructureSplitting COMMAND test_structure_splitting)

add_executable(test_segmentation test_segmentation.cpp effedup.cpp)
target_compile_definitions(test_segmentation PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_segmentation PRIVATE ${LLKA_LIB_LINK})
add_test(NAME Segmentation COMMAND test_segmentation)

add_executable(test_extract_backbone test_extract_backbone.cpp effedup.cpp)
target_compile_definitions(test_extract_backbone PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_extract_backbone ${LLKA_LIB_LINK})
add_test(NAME ExtractBackbone COMMAND test_extract_backbone)

add_executable(test_nucleotide test_nucleotide.cpp effedup.cpp)
target_compile_definitions(test_nucleotide PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_nucleotide ${LLKA_LIB_LINK})
add_test(NAME Nucleotide COMMAND test_nucleotide)

add_executable(test_ntc test_ntc.cpp effedup.cpp)
target_compile_definitions(test_ntc PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_ntc ${LLKA_LIB_LINK})
add_test(NAME NtC COMMAND test_ntc)

add_executable(test_superposition test_superposition.cpp effedup.cpp)
target_compile_definitions(test_superposition PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_superposition ${LLKA_LIB_LINK})
add_test(NAME Superposition COMMAND test_superposition)

add_executable(test_connectivity test_connectivity.cpp effedup.cpp)
target_compile_definitions(test_connectivity PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_connectivity ${LLKA_LIB_LINK})
add_test(NAME Connectivity COMMAND test_connectivity)

add_executable(test_similarity test_similarity.cpp effedup.cpp)
target_compile_definitions(test_similarity PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_similarity ${LLKA_LIB_LINK})
add_test(NAME Similarity COMMAND test_similarity)

add_executable(test_measurements test_measurements.cpp effedup.cpp)
target_compile_definitions(test_measurements PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_measurements ${LLKA_LIB_LINK})
add_test(NAME Measurements COMMAND test_measurements)

if ((NOT EMSCRIPTEN) OR (EMX_JS_BUILD_MODE STREQUAL "Node"))
    add_executable(test_classification test_classification.cpp effedup.cpp)
    target_compile_definitions(test_classification PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
    target_link_libraries(test_classification ${LLKA_LIB_LINK})
    add_test(NAME Classification COMMAND test_classification)
    # Get definition files necessary to initialize classification context
    add_custom_command(
        TARGET test_classification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/clusters.csv"
        "${CMAKE_CURRENT_BINARY_DIR}/clusters.csv"
    )
    add_custom_command(
        TARGET test_classification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/confals.csv"
        "${CMAKE_CURRENT_BINARY_DIR}/confals.csv"
    )
    add_custom_command(
        TARGET test_classification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/golden_steps.csv"
        "${CMAKE_CURRENT_BINARY_DIR}/golden_steps.csv"
    )
    add_custom_command(
        TARGET test_classification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/nu_angles.csv"
        "${CMAKE_CURRENT_BINARY_DIR}/nu_angles.csv"
    )
    add_custom_command(
        TARGET test_classification POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/confal_percentiles.csv"
        "${CMAKE_CURRENT_BINARY_DIR}/confal_percentiles.csv"
    )


    add_executable(test_minicif test_minicif.cpp effedup.cpp)
    target_compile_definitions(test_minicif PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
    # Make sure that the test does not run out of memory when ran with Node.js
    if (EMSCRIPTEN)
        set_target_properties(
            test_minicif
            PROPERTIES
                LINK_FLAGS "-sALLOW_MEMORY_GROWTH -fexceptions"
        )
        target_link_libraries(test_minicif PRIVATE ${LLKA_LIB_LINK} nodefs.js noderawfs.js)
    else()
        target_link_libraries(test_minicif PRIVATE ${LLKA_LIB_LINK})
    endif ()
    add_test(NAME MiniCif COMMAND test_minicif)
    # Get example CIF files
    add_custom_command(
        TARGET test_minicif POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/test_cifs/1BNA.cif"
        "${CMAKE_CURRENT_BINARY_DIR}/1BNA.cif"
    )
    add_custom_command(
        TARGET test_minicif POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/test_cifs/1BNA_broken.cif"
        "${CMAKE_CURRENT_BINARY_DIR}/1BNA_broken.cif"
    )
    add_custom_command(
        TARGET test_minicif POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/test_cifs/quote_in_quotes.cif"
        "${CMAKE_CURRENT_BINARY_DIR}/quote_in_quotes.cif"
    )
    add_custom_command(
        TARGET test_minicif POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../assets/test_cifs/unterminated_quote.cif"
        "${CMAKE_CURRENT_BINARY_DIR}/unterminated_quote.cif"
    )
endif ()

add_executable(test_cpp_interface test_cpp_interface.cpp effedup.cpp)
target_compile_definitions(test_cpp_interface PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(test_cpp_interface ${LLKA_LIB_LINK} ${CPP_EXTRA_LINK})
add_test(NAME CppInterface COMMAND test_cpp_interface)

add_executable(test_cpp_interface_2 test_cpp_interface_2.cpp effedup.cpp)
target_compile_definitions(test_cpp_interface_2 PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
add_test(NAME CppInterface2 COMMAND test_cpp_interface_2)
if (EMSCRIPTEN)
    set_target_properties(
        test_cpp_interface_2
        PROPERTIES
            LINK_FLAGS "-sALLOW_MEMORY_GROWTH -fexceptions"
    )
    target_link_libraries(test_cpp_interface_2 PRIVATE ${LLKA_LIB_LINK} ${CPP_EXTRA_LINK} nodefs.js noderawfs.js)
else ()
    target_link_libraries(test_cpp_interface_2 PRIVATE ${LLKA_LIB_LINK} ${CPP_EXTRA_LINK})
endif ()

add_executable(test_cpp_interface_3 test_cpp_interface_3.cpp effedup.cpp)
target_compile_definitions(test_cpp_interface_3 PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
add_test(NAME CppInterface3 COMMAND test_cpp_interface_3)
target_link_libraries(test_cpp_interface_3 PRIVATE ${LLKA_LIB_LINK} ${CPP_EXTRA_LINK})
