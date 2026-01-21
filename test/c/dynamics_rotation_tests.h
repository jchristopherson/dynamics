#ifndef DYNAMICS_ROTATION_TESTS_H_
#define DYNAMICS_ROTATION_TESTS_H_

#include <stdbool.h>

bool c_test_rotation_x();
bool c_test_rotation_y();
bool c_test_rotation_z();
bool c_test_rotation();
bool c_test_rotation();
bool c_test_acceleration_transform();
bool c_test_velocity_transform();
bool c_test_quaternion_from_array();
bool c_test_quaternion_from_angle_axis();
bool c_test_quaternion_from_matrix();
bool c_test_quaternion_normalize();
bool c_test_quaternion_add();
bool c_test_quaternion_subtract();
bool c_test_quaternion_multiply();
bool c_test_quaternion_divide();
bool c_test_quaternion_abs();
bool c_test_quaternion_to_matrix();

#endif