include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)

# FRF Example 1
add_executable(
    frf_example_1
    frf_example_1.f90
    example_models.f90
)
link_library(frf_example_1 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(frf_example_1 ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# FRF Sweep Example 1
add_executable(
    frf_sweep_example_1
    frf_sweep_example_1.f90
)
link_library(frf_sweep_example_1 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(frf_sweep_example_1 ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
