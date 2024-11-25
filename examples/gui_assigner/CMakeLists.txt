if (BUILD_STATIC_LIBRARY)
    set(LLKA_LIB_LINK libLLKA_STATIC)
else ()
    set(LLKA_LIB_LINK libLLKA_SHARED)
endif ()

find_package(Qt6 COMPONENTS Core Gui Widgets)
qt_standard_project_setup()

set(
    gui_assigner_SRCS
    src/assignment.cpp
    src/main.cpp
    src/worker.cpp
    src/ui/assignment_in_progress_dlg.cpp
    src/ui/assignment_parameters_dlg.cpp
    src/ui/detailed_step_info.cpp
    src/ui/detailed_step_info_dlg.cpp
    src/ui/main_window.cpp
    src/ui/setup_tracepoints_dlg.cpp
    src/ui/view_last_trace_dlg.cpp
    src/ui/models/ntc_assignment_table_model.cpp
)

if (WIN32)
    add_executable(gui_assigner WIN32 ${gui_assigner_SRCS})
else ()
    add_executable(gui_assigner ${gui_assigner_SRCS})
endif ()
target_compile_definitions(gui_assigner PRIVATE ${LIBLLKA_GLOBAL_DEFINITIONS} ${LIBLLKA_PLATFORM_DEFINITIONS})
target_link_libraries(
    gui_assigner
    PRIVATE ${LLKA_LIB_LINK}
    PRIVATE Qt6::Core
    PRIVATE Qt6::Gui
    PRIVATE Qt6::Widgets
)
