# Get Dependencies
include(FetchContent)
find_package(fplot QUIET)
find_package(diffeq QUIET)

if (NOT fplot_FOUND)
    FetchContent_Declare(
        fplot
        GIT_TAG "origin/master"
        GIT_REPOSITORY https://github.com/jchristopherson/fplot
    )
    FetchContent_MakeAvailable(fplot)
    set(FPLOT_LIBRARIES fplot)
else()
    set(FPLOT_LIBRARIES fplot::fplot)
endif()

if (NOT diffeq_FOUND)
    message(STATUS "DIFFEQ could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        diffeq
        GIT_REPOSITORY "https://github.com/jchristopherson/diffeq"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(diffeq)
    set(DIFFEQ_LIBRARIES diffeq)
else()
    set(DIFFEQ_LIBRARIES diffeq::diffeq)
endif()

# FRF Sweep Example 1
add_executable(
    frf_sweep_example_1
    frf_sweep_example_1.f90
)
target_link_libraries(
    frf_sweep_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# FRF Proportional Damping Example 1
add_executable(
    frf_proportional_example_1
    frf_proportional_example_1.f90
)
target_link_libraries(
    frf_proportional_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Rotation Example 1
add_executable(
    rotation_example_1
    rotation_example_1.f90
)
target_link_libraries(
    rotation_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Kinematics Example 1
add_executable(
    kinematics_example_1
    kinematics_example_1.f90
)
target_link_libraries(kinematics_example_1 ${PROJECT_NAME})

# Linkage Example 1
add_executable(
    linkage_example_1
    linkage_example_1.f90
)
target_link_libraries(linkage_example_1 ${PROJECT_NAME} ${COLLECTIONS_LIBRARIES})

# SISO FRF Example 1
add_executable(
    siso_frf_example_1
    siso_frf_example_1.f90
)
target_link_libraries(
    siso_frf_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
)

# MIMO FRF Example 1
add_executable(
    mimo_frf_example_1
    mimo_frf_example_1.f90
)
target_link_libraries(
    mimo_frf_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
)

# FRF Fitting Example 1
add_executable(
    frf_fit_example_1
    frf_fit_example_1.f90
)
target_link_libraries(
    frf_fit_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Free Vibration Example
add_executable(
    free_vibration_example
    free_vibration_example.f90
)
target_link_libraries(
    free_vibration_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Jacobian Example 1
add_executable(
    jacobian_example_1
    jacobian_example_1.f90
)
target_link_libraries(
    jacobian_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# SISO Fitting Example
add_executable(
    siso_lsq_fit_example
    siso_lsq_fit_example.f90
)
target_link_libraries(
    siso_lsq_fit_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
)

# Poincare Example
add_executable(
    poincare_example
    poincare_example.f90
)
target_link_libraries(
    poincare_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
)

# C DH_Table Example
if (${BUILD_DYNAMICS_C_INTERFACE})
    add_executable(
        c_dh_table_example_1
        c_dh_table_example_1.c
    )
    target_link_libraries(
        c_dh_table_example_1
        ${PROJECT_NAME}
    )
    target_include_directories(
        c_dh_table_example_1
        PUBLIC
        ${CMAKE_SOURCE_DIR}/src/c
    )
endif()

# PID Control Example 1
add_executable(
    pid_control_example_1
    pid_control_example_1.f90
)
target_link_libraries(
    pid_control_example_1
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# C LTI Example 1
if (${BUILD_DYNAMICS_C_INTERFACE})
    add_executable(
        c_lti_example_1
        c_lti_example_1.c
    )
    target_link_libraries(
        c_lti_example_1
        ${PROJECT_NAME}
    )
    target_include_directories(
        c_lti_example_1
        PUBLIC
        ${CMAKE_SOURCE_DIR}/src/c
    )
endif()

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET frf_sweep_example_1
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:frf_sweep_example_1> 
        $<TARGET_FILE_DIR:frf_sweep_example_1>
        COMMAND_EXPAND_LISTS
    )
endif()