# Define the source files
set(dynamics_test_sources
    dynamics_tests.f90
    dynamics_frf_tests.f90
    dynamics_structures_tests.f90
    dynamics_kinematics_tests.f90
    dynamics_vibrations_tests.f90
    dynamics_stability_tests.f90
    dynamics_transfer_function_tests.f90
    dynamics_state_space_tests.f90
    dynamics_system_id_tests.f90
    dynamics_rotation_tests.f90
)

# Dependencies
include(FetchContent)
find_package(fortran_test_helper QUIET)
find_package(ferror QUIET)
find_package(linalg QUIET)
find_package(diffeq QUIET)
find_package(fstats QUIET)
find_package(collections QUIET)
if (NOT fortran_test_helper_FOUND)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
else()
    set(FERROR_LIBRARIES ferror::ferror)
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

# Build the tests
add_executable(dynamics_tests ${dynamics_test_sources})
target_link_libraries(
    dynamics_tests
    dynamics
    ${FTH_LIBRARIES}
    ${FERROR_LIBRARIES}
    ${LINALG_LIBRARIES}
    ${FSTATS_LIBRARIES}
    ${COLLECTIONS_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
)
add_test(
    NAME dynamics_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:dynamics_tests>
)

# Test the C API
if (BUILD_DYNAMICS_C_INTERFACE)
    add_subdirectory(c)
endif()

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET dynamics_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:dynamics_tests> 
        $<TARGET_FILE_DIR:dynamics_tests>
        COMMAND_EXPAND_LISTS
    )
endif()
