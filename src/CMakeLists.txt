# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DYNAMICS_SOURCES
    ${dir}/dynamics.f90
    ${dir}/dynamics_error_handling.f90
    ${dir}/dynamics_frequency_response.f90
    ${dir}/dynamics_rotation.f90
    ${dir}/dynamics_structural.f90
    ${dir}/dynamics_kinematics.f90
    ${dir}/dynamics_vibrations.f90
    ${dir}/dynamics_helper.f90
    ${dir}/dynamics_stability.f90
    ${dir}/dynamics_controls.f90
    ${dir}/dynamics_system_id.f90
)
set(DYNAMICS_SOURCES ${DYNAMICS_SOURCES} PARENT_SCOPE)