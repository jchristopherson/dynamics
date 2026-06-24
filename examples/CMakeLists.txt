# Get Dependencies
include(FetchContent)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
find_package(ferror QUIET)
find_package(fplot QUIET)
find_package(fstats QUIET)
find_package(collections QUIET)
find_package(diffeq QUIET)
find_package(linalg QUIET)
find_package(spectrum QUIET)
find_package(geompack QUIET)
find_package(fstring QUIET)
if (NOT ferror_FOUND)
    FetchContent_Declare(
        ferror
        GIT_TAG "origin/master"
        GIT_REPOSITORY https://github.com/jchristopherson/ferror
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
else()
    set(FERROR_LIBRARIES ferror::ferror)
endif()
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

if (NOT fstats_FOUND)
    message(STATUS "FSTATS could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        fstats
        GIT_REPOSITORY "https://github.com/jchristopherson/fstats"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(fstats)
    set(FSTATS_LIBRARIES fstats)
else()
    set(FSTATS_LIBRARIES fstats::fstats)
endif()

if (NOT collections_FOUND)
    message(STATUS "COLLECTIONS could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        collections
        GIT_REPOSITORY "https://github.com/jchristopherson/collections"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(collections)
    set(COLLECTIONS_LIBRARIES collections)
else()
    set(COLLECTIONS_LIBRARIES collections::collections)
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

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
else()
    set(LINALG_LIBRARIES linalg::linalg)
endif()

if (NOT spectrum_FOUND)
    message(STATUS "SPECTRUM could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        spectrum
        GIT_REPOSITORY "https://github.com/jchristopherson/spectrum"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(spectrum)
    set(SPECTRUM_LIBRARIES spectrum)
else()
    set(SPECTRUM_LIBRARIES spectrum::spectrum)
endif()

if (NOT geompack_FOUND)
    message(STATUS "GEOMPACK could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        geompack
        GIT_REPOSITORY "https://github.com/jchristopherson/geompack"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(geompack)
    set(GEOMPACK_LIBRARIES geompack)
else()
    set(GEOMPACK_LIBRARIES geompack::geompack)
endif()

if (NOT fstring_FOUND)
    message(STATUS "FSTRING could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        fstring
        GIT_REPOSITORY "https://github.com/jchristopherson/fstring"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(fstring)
    set(FSTRING_LIBRARIES fstring)
else()
    set(FSTRING_LIBRARIES fstring::fstring)
endif()

FetchContent_Declare(
    forcolormap
    GIT_REPOSITORY "https://github.com/vmagnin/forcolormap"
    GIT_TAG main
)
FetchContent_MakeAvailable(forcolormap)

# FRF Sweep Example 1
add_executable(
    frf_sweep_example_1
    frf_sweep_example_1.f90
)
target_link_libraries(
    frf_sweep_example_1
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# FRF Proportional Damping Example 1
add_executable(
    frf_proportional_example_1
    frf_proportional_example_1.f90
)
target_link_libraries(
    frf_proportional_example_1
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# Rotation Example 1
add_executable(
    rotation_example_1
    rotation_example_1.f90
)
target_link_libraries(
    rotation_example_1
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
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
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
    ${LINALG_LIBRARIES}
    ${SPECTRUM_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# MIMO FRF Example 1
add_executable(
    mimo_frf_example_1
    mimo_frf_example_1.f90
)
target_link_libraries(
    mimo_frf_example_1
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
    ${LINALG_LIBRARIES}
    ${SPECTRUM_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# FRF Fitting Example 1
add_executable(
    frf_fit_example_1
    frf_fit_example_1.f90
)
target_link_libraries(
    frf_fit_example_1
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# Free Vibration Example
add_executable(
    free_vibration_example
    free_vibration_example.f90
)
target_link_libraries(
    free_vibration_example
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# Jacobian Example 1
add_executable(
    jacobian_example_1
    jacobian_example_1.f90
)
target_link_libraries(
    jacobian_example_1
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)

# SISO Fitting Example
add_executable(
    siso_lsq_fit_example
    siso_lsq_fit_example.f90
)
target_link_libraries(
    siso_lsq_fit_example
    ${PROJECT_NAME}
    ${FERROR_LIBRARIES}
    ${FPLOT_LIBRARIES}
    ${FSTATS_LIBRARIES}
    ${COLLECTIONS_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
    ${LINALG_LIBRARIES}
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
)
if (OpenMP_Fortran_FOUND)
    target_link_libraries(siso_lsq_fit_example OpenMP::OpenMP_Fortran)
endif()

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
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
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
    ${GEOMPACK_LIBRARIES}
    ${FSTRING_LIBRARIES}
    $<BUILD_INTERFACE:forcolormap>
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

if (OpenMP_Fortran_FOUND)
    target_link_libraries(frf_sweep_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(frf_proportional_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(rotation_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(kinematics_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(linkage_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(siso_frf_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(mimo_frf_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(frf_fit_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(free_vibration_example OpenMP::OpenMP_Fortran)
    target_link_libraries(jacobian_example_1 OpenMP::OpenMP_Fortran)
    target_link_libraries(siso_lsq_fit_example OpenMP::OpenMP_Fortran)
    target_link_libraries(poincare_example OpenMP::OpenMP_Fortran)
    target_link_libraries(pid_control_example_1 OpenMP::OpenMP_Fortran)
    
    if (${BUILD_DYNAMICS_C_INTERFACE})
        target_link_libraries(c_dh_table_example_1 OpenMP::OpenMP_Fortran)
        target_link_libraries(c_lti_example_1 OpenMP::OpenMP_Fortran)
    endif()
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