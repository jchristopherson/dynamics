set(dir ${CMAKE_CURRENT_LIST_DIR})

set(DYNAMICS_C_SOURCES
    ${dir}/dynamics_c_api.f90
    ${dir}/dynamics.c
)
set(DYNAMICS_C_SOURCES ${DYNAMICS_C_SOURCES} PARENT_SCOPE)