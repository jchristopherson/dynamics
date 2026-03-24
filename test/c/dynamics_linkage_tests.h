#ifndef DYNAMICS_LINKAGE_TESTS_H_
#define DYNAMICS_LINKAGE_TESTS_H_

#include <stdbool.h>
#include "dynamics.h"

c_serial_linkage build_serial_linkage_1();
bool c_test_serial_linkage_forward_kinematics();
bool c_test_serial_linkage_jacobian();
bool c_test_serial_linkage_inverse_kinematics();

#endif
