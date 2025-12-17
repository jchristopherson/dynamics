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

# C API
add_subdirectory(c)
if (${BUILD_DYNAMICS_C_INTERFACE})
    list(APPEND DYNAMICS_SOURCES ${DYNAMICS_C_SOURCES})
endif()

# Dependencies
include(FetchContent)

find_package(ferror QUIET)
find_package(linalg QUIET)
find_package(nonlin QUIET)
find_package(fstats QUIET)
find_package(diffeq QUIET)
find_package(peaks QUIET)
find_package(spectrum QUIET)
find_package(collections QUIET)

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
)
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    ${LINALG_LIBRARIES}
    ${FSTATS_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
    ${PEAKS_LIBRARIES}
    ${SPECTRUM_LIBRARIES}
    ${COLLECTIONS_LIBRARIES}
)

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

if (${BUILD_FERROR} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FERROR_LIBRARIES})
endif()

if (${BUILD_LINALG} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${LINALG_LIBRARIES})
endif()

if (${BUILD_NONLIN} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${NONLIN_LIBRARIES})
endif()

if (${BUILD_FSTATS} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FSTATS_LIBRARIES})
endif()

if (${BUILD_DIFFEQ} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${DIFFEQ_LIBRARIES})
endif()

if (${BUILD_PEAKS} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${PEAKS_LIBRARIES})
endif()

if (${BUILD_SPECTRUM} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${SPECTRUM_LIBRARIES})
endif()

if (${BUILD_COLLECTIONS} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${COLLECTIONS_LIBRARIES})
endif()
