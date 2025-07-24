#ifndef DYNAMICS_VIBRATIONS_TESTS_H_
#define DYNAMICS_VIBRATIONS_TESTS_H_

#include <stdbool.h>

bool c_test_q_factor();
bool c_test_bandwidth();
bool c_test_log_decrement();
bool c_test_damping_from_decrement();
bool c_test_find_free_rsp_props();
bool c_test_rise_time();
bool c_test_step_response();
bool c_test_settling_amplitude();
bool c_test_damping_from_overshoot();

#endif