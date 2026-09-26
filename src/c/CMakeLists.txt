set(dir ${CMAKE_CURRENT_LIST_DIR})

set(DYNAMICS_C_SOURCES
    ${dir}/dynamics_c_api.f90
    ${dir}/dynamics_c_types.f90
    ${dir}/dynamics_c_support.f90
    ${dir}/dynamics_c_rotation.f90
    ${dir}/dynamics_c_stability.f90
    ${dir}/dynamics_c_kinematics.f90
    ${dir}/dynamics_c_frequency.f90
    ${dir}/dynamics_c_system_id.f90
    ${dir}/dynamics_c_quaternions.f90
    ${dir}/dynamics_c_geometry.f90
    ${dir}/dynamics_c_maps.f90
    ${dir}/dynamics_c_linkage.f90
    ${dir}/dynamics_c_parallel_linkage.f90
    ${dir}/dynamics_c_controls.f90
    ${dir}/dynamics_c_structural.f90
    ${dir}/dynamics_c_vibrations.f90
    ${dir}/dynamics_c_variational.f90
    ${dir}/dynamics.c
)
set(DYNAMICS_C_SOURCES ${DYNAMICS_C_SOURCES} PARENT_SCOPE)