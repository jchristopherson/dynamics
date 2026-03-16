#ifndef DYNAMICS_KINEMATICS_TESTS_H_
#define DYNAMICS_KINEMATICS_TESTS_H_

#include <stdbool.h>

bool c_test_forward_kinematics();
bool c_test_inverse_kinematics();

bool c_test_define_link_csys();
bool c_test_define_csys();
bool c_test_build_dh_table();

#endif
