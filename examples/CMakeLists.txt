include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)

# FRF Sweep Example 1
add_executable(
    frf_sweep_example_1
    frf_sweep_example_1.f90
)
link_library(frf_sweep_example_1 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(frf_sweep_example_1 ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# FRF Proportional Damping Example 1
add_executable(
    frf_proportional_example_1
    frf_proportional_example_1.f90
)
link_library(frf_proportional_example_1 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(frf_proportional_example_1 ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Rotation Example 1
add_executable(
    rotation_example_1
    rotation_example_1.f90
)
link_library(rotation_example_1 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rotation_example_1 ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Kinematics Example 1
add_executable(
    kinematics_example_1
    kinematics_example_1.f90
)
link_library(kinematics_example_1 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
