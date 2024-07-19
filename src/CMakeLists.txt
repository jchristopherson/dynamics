# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DYNAMICS_SOURCES
    ${dir}/dynamics.f90
    ${dir}/dynamics_error_handling.f90
    ${dir}/dynamics_frequency_response.f90
    ${dir}/dynamics_rotation.f90
    ${dir}/dynamics_structural.f90
)
set(DYNAMICS_SOURCES ${DYNAMICS_SOURCES} PARENT_SCOPE)