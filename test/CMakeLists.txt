# Define the source files
set(dynamics_test_sources
    dynamics_tests.f90
    dynamics_frf_tests.f90
    dynamics_structures_tests.f90
    dynamics_kinematics_tests.f90
    dynamics_vibrations_tests.f90
    dynamics_stability_tests.f90
    dynamics_transfer_function_tests.f90
    dynamics_state_space_tests.f90
    dynamics_system_id_tests.f90
    dynamics_rotation_tests.f90
    dynamics_geometry_tests.f90
    dynamics_linkage_tests.f90
    dynamics_parallel_linkage_tests.f90
)

# Dependencies
find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    set(_linalg_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(fortran_test_helper)
    set(BUILD_TESTING ${_linalg_build_testing})
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

# Build the tests
add_executable(dynamics_tests ${dynamics_test_sources})
target_link_libraries(
    dynamics_tests
    dynamics
    ${FTH_LIBRARIES}
)
add_test(
    NAME dynamics_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:dynamics_tests>
)

# Test the C API
if (BUILD_DYNAMICS_C_INTERFACE)
    add_subdirectory(c)
endif()

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET dynamics_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:dynamics_tests> 
        $<TARGET_FILE_DIR:dynamics_tests>
        COMMAND_EXPAND_LISTS
    )
endif()
