# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    dynamics
    LANGUAGES Fortran C CXX
    VERSION 1.2.0
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED TRUE)

# Confgiure everything
add_subdirectory(configure)

# Build the C API
option(BUILD_DYNAMICS_C_INTERFACE "Build DYNAMICS C interface")
message(STATUS "Build DYNAMICS C interface: ${BUILD_DYNAMICS_C_INTERFACE}")

# Source
add_subdirectory(src)

# Installation Instructions
install(
  EXPORT ${PROJECT_NAME}-targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build DYNAMICS tests: ${BUILD_TESTING}")
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

# Additional Compiler Options
if (${BUILD_DYNAMICS_C_INTERFACE} AND ${BUILD_SHARED_LIBS})
    target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)
endif()
