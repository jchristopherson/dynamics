# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Define Source Files
set(dynamics_c_sources
    dynamics_c_tests.c
    dynamics_c_test_helper.c
    dynamics_vibrations_tests.c
    dynamics_rotation_tests.c
    dynamics_stability_tests.c
)

# Build the tests
add_executable(dynamics_c_tests ${dynamics_c_sources})
link_library(dynamics_c_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
target_include_directories(dynamics_c_tests PRIVATE
    "${PROJECT_SOURCE_DIR}/interfaces/c"
)
target_include_directories(dynamics_c_tests PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
)
add_test(
    NAME dynamics_c_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:dynamics_c_tests>
)