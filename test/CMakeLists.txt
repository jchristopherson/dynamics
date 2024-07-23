# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Dependencies
add_subdirectory(fortran_test_helper)

# Define the source files
set(dynamics_test_sources
    dynamics_tests.f90
    dynamics_frf_tests.f90
    dynamics_structures_tests.f90
)

# Build the tests
add_executable(dynamics_tests ${dynamics_test_sources})
link_library(dynamics_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(dynamics_tests ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME dynamics_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:dynamics_tests>
)
