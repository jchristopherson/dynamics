#ifndef DYNAMICS_FRF_TESTS_H_
#define DYNAMICS_FRF_TESTS_H_

#include <stdbool.h>

bool c_test_frequency_response();
bool c_test_general_damping_frf();
bool c_test_dynamic_stiffness_dense();
bool c_test_modal_response();
bool c_test_frf_sweep();
bool c_test_frf_fit();
bool c_test_siso_frf();
bool c_test_siso_lsq_fit();

#endif