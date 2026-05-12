#include "dynamics_vibrations_tests.h"
#include "dynamics_rotation_tests.h"
#include "dynamics_stability_tests.h"
#include "dynamics_kinematics_tests.h"
#include "dynamics_frf_tests.h"
#include "dynamics_geometry_tests.h"
#include "dynamics_linkage_tests.h"
#include "dynamics_transfer_function_tests.h"
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
    if (!c_test_quaternion_exp()) flag = 36;
    if (!c_test_quaternion_log()) flag = 37;
    if (!c_test_quaternion_pow()) flag = 38;
    if (!c_test_quaternion_dot_product()) flag = 39;
    if (!c_test_quaternion_roll_pitch_yaw()) flag = 40;

    if (!c_test_line_from_2_points()) flag = 41;
    if (!c_test_line_from_2_planes()) flag = 42;
    if (!c_test_line_eval()) flag = 43;
    if (!c_test_plane_from_3_points()) flag = 44;
    if (!c_test_plane_from_point_and_normal()) flag = 45;
    if (!c_test_plane_normal()) flag = 46;
    if (!c_test_is_parallel_vectors()) flag = 47;
    if (!c_test_is_parallel_planes()) flag = 48;
    if (!c_test_is_point_on_plane()) flag = 49;
    if (!c_test_point_to_line_distance()) flag = 50;
    if (!c_test_is_point_on_line()) flag = 51;
    if (!c_test_point_to_plane_distance()) flag = 52;
    if (!c_test_vector_plane_projection()) flag = 53;
    if (!c_test_fit_line_to_many_points()) flag = 54;
    if (!c_test_flip_plane_normal()) flag = 55;

    if (!c_test_plucker_line_from_2_points()) flag = 56;
    if (!c_test_plucker_line_from_line()) flag = 57;
    if (!c_test_plucker_line_from_2_planes()) flag = 58;
    if (!c_test_plucker_line_matmul()) flag = 59;

    if (!c_test_define_link_csys()) flag = 60;
    if (!c_test_define_link_csys()) flag = 61;
    if (!c_test_build_dh_table()) flag = 62;
    
    if (!c_test_line_common_normal()) flag = 63;
    if (!c_test_do_lines_intersect()) flag = 64;

    if (!c_test_serial_linkage_forward_kinematics()) flag = 65;
    if (!c_test_serial_linkage_jacobian()) flag = 66;
    if (!c_test_serial_linkage_inverse_kinematics()) flag = 67;

    if (!c_test_tf_evaluate()) flag = 68;
    if (!c_test_ccf_form()) flag = 69;
    if (!c_test_ocf_form()) flag = 70;
    if (!c_test_poles_zeros()) flag = 71;
    if (!c_test_tf_multiply()) flag = 72;
    if (!c_test_state_space_initialize()) flag = 73;

    return flag;
}
