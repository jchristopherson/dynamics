#include "dynamics_vibrations_tests.h"
#include "dynamics_rotation_tests.h"
#include "dynamics_stability_tests.h"
#include "dynamics_kinematics_tests.h"
#include "dynamics_frf_tests.h"
#include <stdbool.h>

int main()
{
    // Local Variables
    int flag;

    // Initialization
    flag = 0;

    // Tests
    if (!c_test_q_factor()) flag = 1;
    if (!c_test_bandwidth()) flag = 2;
    if (!c_test_log_decrement()) flag = 3;
    if (!c_test_damping_from_decrement()) flag = 4;
    if (!c_test_find_free_rsp_props()) flag = 5;
    if (!c_test_rise_time()) flag = 6;
    if (!c_test_step_response()) flag = 7;
    if (!c_test_settling_amplitude()) flag = 8;
    if (!c_test_damping_from_overshoot()) flag = 9;

    if (!c_test_rotation_x()) flag = 10;
    if (!c_test_rotation_y()) flag = 11;
    if (!c_test_rotation_z()) flag = 12;
    if (!c_test_rotation()) flag = 13;
    if (!c_test_acceleration_transform()) flag = 14;
    if (!c_test_velocity_transform()) flag = 15;

    if (!c_test_determine_local_stability()) flag = 16;

    if (!c_test_forward_kinematics()) flag = 17;
    if (!c_test_inverse_kinematics()) flag = 18;

    if (!c_test_frequency_response()) flag = 19;
    if (!c_test_modal_response()) flag = 20;
    if (!c_test_frf_sweep()) flag = 21;
    if (!c_test_frf_fit()) flag = 22;
    if (!c_test_siso_frf()) flag = 23;
    if (!c_test_siso_lsq_fit()) flag = 24;

    if (!c_test_quaternion_from_array()) flag = 25;
    if (!c_test_quaternion_from_angle_axis()) flag = 26;
    if (!c_test_quaternion_from_matrix()) flag = 27;
    if (!c_test_quaternion_normalize()) flag = 28;
    if (!c_test_quaternion_add()) flag = 29;
    if (!c_test_quaternion_subtract()) flag = 30;
    if (!c_test_quaternion_multiply()) flag = 31;
    if (!c_test_quaternion_divide()) flag = 32;
    if (!c_test_quaternion_abs()) flag = 33;
    if (!c_test_quaternion_to_matrix()) flag = 34;
    if (!c_test_to_angle_axis()) flag = 35;

    return flag;
}
