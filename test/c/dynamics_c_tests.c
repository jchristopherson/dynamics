#include "dynamics_vibrations_tests.h"
#include "dynamics_rotation_tests.h"
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

    return flag;
}
