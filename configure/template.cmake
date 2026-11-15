@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(OpenMP COMPONENTS Fortran)
find_dependency(linalg 3.0 QUIET)
find_dependency(collections 1.2 QUIET)
find_dependency(diffeq 3.0 QUIET)
find_dependency(fstats 2.0 QUIET)
find_dependency(nonlin 3.0 QUIET)
find_dependency(peaks QUIET)
find_dependency(spectrum 1.5 QUIET)
find_dependency(autodiff 1.0 QUIET)

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()