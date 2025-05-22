# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    dynamics
    LANGUAGES Fortran
    VERSION 1.0.11
)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
add_subdirectory(dependencies)

# Source
add_subdirectory(src)
add_fortran_library(
    ${PROJECT_NAME}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_INCLUDEDIR}
    ${PROJECT_VERSION}
    ${PROJECT_VERSION_MAJOR}
    ${DYNAMICS_SOURCES}
)
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${linalg_LIBRARY} ${linalg_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${spectrum_LIBRARY} ${spectrum_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${diffeq_LIBRARY} ${diffeq_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${nonlin_LIBRARY} ${nonlin_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${fstats_LIBRARY} ${fstats_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${peaks_LIBRARY} ${peaks_INCLUDE_DIR})

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# Examples
option(BUILD_DYNAMICS_EXAMPLES "Build DYNAMICS examples")
message(STATUS "Build DYNAMICS examples: ${BUILD_DYNAMICS_EXAMPLES}")
if (BUILD_DYNAMICS_EXAMPLES)
    add_subdirectory(examples)
endif()