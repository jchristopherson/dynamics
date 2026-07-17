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
    ${dir}/dynamics_quaternions.f90
    ${dir}/dynamics_maps.f90
    ${dir}/dynamics_geometry.f90
    ${dir}/dynamics_rigid_bodies.f90
    ${dir}/dynamics_linkage.f90
)

# C API
add_subdirectory(c)
if (${BUILD_DYNAMICS_C_INTERFACE})
    list(APPEND DYNAMICS_SOURCES ${DYNAMICS_C_SOURCES})
endif()

# Dependencies
include(FetchContent)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
find_package(ferror QUIET)
find_package(linalg 1.8.10 QUIET)
find_package(nonlin 2.1.4 QUIET)
find_package(fstats 1.5.6 QUIET)
find_package(diffeq 1.1.5 QUIET)
find_package(peaks QUIET)
find_package(spectrum QUIET)
find_package(collections 1.0.9 QUIET)

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
    set(BUILD_FERROR TRUE)
else()
    set(FERROR_LIBRARIES ferror::ferror)
    set(BUILD_FERROR FALSE)
endif()

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
    set(BUILD_LINALG TRUE)
else()
    set(LINALG_LIBRARIES linalg::linalg)
    set(BUILD_LINALG FALSE)
endif()

if (NOT nonlin_FOUND)
    message(STATUS "NONLIN could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        nonlin
        GIT_REPOSITORY "https://github.com/jchristopherson/nonlin"
    )
    FetchContent_MakeAvailable(nonlin)
    set(NONLIN_LIBRARIES nonlin)
    set(BUILD_NONLIN TRUE)
else()
    set(NONLIN_LIBRARIES nonlin::nonlin)
    set(BUILD_NONLIN FALSE)
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
    set(BUILD_FSTATS TRUE)
else()
    set(FSTATS_LIBRARIES fstats::fstats)
    set(BUILD_FSTATS FALSE)
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
    set(BUILD_DIFFEQ TRUE)
else()
    set(DIFFEQ_LIBRARIES diffeq::diffeq)
    set(BUILD_DIFFEQ FALSE)
endif()

if (NOT peaks_FOUND)
    message(STATUS "PEAKS could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        peaks
        GIT_REPOSITORY "https://github.com/jchristopherson/peaks"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(peaks)
    set(PEAKS_LIBRARIES peaks)
    set(BUILD_PEAKS TRUE)
else()
    set(PEAKS_LIBRARIES peaks::peaks)
    set(BUILD_PEAKS FALSE)
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
    set(BUILD_SPECTRUM TRUE)
else()
    set(SPECTRUM_LIBRARIES spectrum::spectrum)
    set(BUILD_SPECTRUM FALSE)
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
    set(BUILD_COLLECTIONS TRUE)
else()
    set(COLLECTIONS_LIBRARIES collections::collections)
    set(BUILD_COLLECTIONS FALSE)
endif()

# Build
add_library(${PROJECT_NAME} ${DYNAMICS_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${FERROR_LIBRARIES}
    ${NONLIN_LIBRARIES}
    ${FSTATS_LIBRARIES}
    ${LINALG_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
    ${PEAKS_LIBRARIES}
    ${SPECTRUM_LIBRARIES}
    ${COLLECTIONS_LIBRARIES}
)
if (OpenMP_Fortran_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USEOPENMP=1)

    if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fopenmp>)
        target_link_options(${PROJECT_NAME} PUBLIC $<$<LINK_LANGUAGE:Fortran>:-fopenmp>)
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
        target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-qopenmp>)
        target_link_options(${PROJECT_NAME} PUBLIC $<$<LINK_LANGUAGE:Fortran>:-qopenmp>)
    elseif (MSVC)
        target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:/openmp>)
        target_link_options(${PROJECT_NAME} PUBLIC $<$<LINK_LANGUAGE:Fortran>:/openmp>)
    endif()
endif()

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(DYNAMICS_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${DYNAMICS_MOD_DIR}")
    file(MAKE_DIRECTORY "${DYNAMICS_MOD_DIR}")
endif()

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${DYNAMICS_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${DYNAMICS_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

# Install
install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${DYNAMICS_MOD_DIR} DESTINATION ${CMAKE_INSTALL_MODULEDIR})

if (${BUILD_SHARED_LIBS})
    if (${BUILD_FERROR})
        install(TARGETS ${FERROR_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_COLLECTIONS})
        install(TARGETS ${COLLECTIONS_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_FSTATS})
        install(TARGETS ${FSTATS_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_DIFFEQ})
        install(TARGETS ${DIFFEQ_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_PEAKS})
        install(TARGETS ${PEAKS_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_NONLIN})
        install(TARGETS ${NONLIN_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_SPECTRUM})
        install(TARGETS ${SPECTRUM_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
else()
    if (${BUILD_FERROR})
        install(TARGETS ${FERROR_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_COLLECTIONS})
        install(TARGETS ${COLLECTIONS_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_FSTATS})
        install(TARGETS ${FSTATS_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_DIFFEQ})
        install(TARGETS ${DIFFEQ_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_PEAKS})
        install(TARGETS ${PEAKS_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_NONLIN})
        install(TARGETS ${NONLIN_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_SPECTRUM})
        install(TARGETS ${SPECTRUM_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
endif()