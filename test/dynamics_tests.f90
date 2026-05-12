program main
    use iso_fortran_env, only : int32
    use dynamics_frf_tests
    use dynamics_structures_tests
    use dynamics_kinematics_tests
    use dynamics_vibrations_tests
    use dynamics_stability_tests
    use dynamics_transfer_function_tests
    use dynamics_state_space_tests
    use dynamics_system_id_tests
    use dynamics_rotation_tests
    use dynamics_geometry_tests
    use dynamics_linkage_tests
    implicit none

    ! Variables
    logical :: check
    integer(int32) :: flag

    ! Initialization
    flag = 0

    ! Tests
    check = test_frf_sweep()
    if (.not.check) flag = 1

    check = test_proportional_damping_frf()
    if (.not.check) flag = 2

    check = test_modal_response()
    if (.not.check) flag = 3

    check = test_beam2d_shape_functions()
    if (.not.check) flag = 4

    check = test_shape_function_derivatives()
    if (.not.check) flag = 5

    check = test_beam2d_strain_displacement()
    if (.not.check) flag = 6

    check = test_beam2d_stiffness_matrix()
    if (.not.check) flag = 7

    check = test_beam2d_mass_matrix()
    if (.not.check) flag = 8

    check = test_beam2d_ext_force()
    if (.not.check) flag = 9

    check = test_boundary_conditions()
    if (.not.check) flag = 10

    check = test_connectivity_matrix()
    if (.not.check) flag = 11

    check = test_forward_kinematics()
    if (.not.check) flag = 12

    check = test_inverse_kinematics()
    if (.not.check) flag = 13

    check = test_boundary_conditions_2()
    if (.not.check) flag = 14

    check = test_beam3d_shape_function_matrix()
    if (.not.check) flag = 15

    check = test_beam3d_strain_displacement()
    if (.not.check) flag = 16

    check = test_beam3d_stiffness_matrix()
    if (.not.check) flag = 17

    check = test_beam3d_mass_matrix()
    if (.not.check) flag = 18

    check = test_frf_fit()
    if (.not.check) flag = 19

    check = test_q_factor()
    if (.not.check) flag = 20

    check = test_bandwidth()
    if (.not.check) flag = 21

    check = test_log_decrement()
    if (.not.check) flag = 22

    check = test_damping_from_decrement()
    if (.not.check) flag = 23

    check = test_find_free_rsp_props()
    if (.not.check) flag = 24

    check = test_step_response()
    if (.not.check) flag = 25

    check = test_settling_amplitude()
    if (.not.check) flag = 26

    check = test_damping_from_overshoot()
    if (.not.check) flag = 27

    check = test_jacobian()
    if (.not.check) flag = 28

    check = test_velocity_matrix()
    if (.not.check) flag = 29

    check = test_acceleration_matrix()
    if (.not.check) flag = 30

    check = test_skew_symmetric()
    if (.not.check) flag = 31

    check = test_local_stability()
    if (.not.check) flag = 32

    check = test_tf_evaluate()
    if (.not.check) flag = 33

    check = test_tf_multiply()
    if (.not.check) flag = 34

    check = test_tf_poles_zeros()
    if (.not.check) flag = 35

    check = test_ccf_form_conversion()
    if (.not.check) flag = 36

    check = test_ocf_form_conversion()
    if (.not.check) flag = 37

    check = test_lti_solve()
    if (.not.check) flag = 38

    check = test_siso_model_fit_least_squares()
    if (.not.check) flag = 39

    check = test_siso_model_fit_least_squares_constrained()
    if (.not.check) flag = 40

    check = test_siso_model_fit_least_squares_multi()
    if (.not.check) flag = 41

    check = test_quaternion_init_1()
    if (.not.check) flag = 42

    check = test_quaternion_init_2()
    if (.not.check) flag = 43

    check = test_quaternion_init_3()
    if (.not.check) flag = 44

    check = test_quaternion_abs()
    if (.not.check) flag = 45

    check = test_quaternion_conjg()
    if (.not.check) flag = 46

    check = test_quaternion_real()
    if (.not.check) flag = 47

    check = test_quaternion_aimag()
    if (.not.check) flag = 48

    check = test_quaternion_multiply()
    if (.not.check) flag = 49

    check = test_quaternion_add()
    if (.not.check) flag = 50

    check = test_quaternion_subtract()
    if (.not.check) flag = 51

    check = test_quaternion_division()
    if (.not.check) flag = 52

    check = test_quaternion_to_matrix()
    if (.not.check) flag = 53

    check = test_quaternion_normalize()
    if (.not.check) flag = 54

    check = test_quaternion_to_array()
    if (.not.check) flag = 55

    check = test_quaternion_inverse()
    if (.not.check) flag = 56

    check = test_quaternion_to_angle_axis()
    if (.not.check) flag = 57

    check = test_quaternion_exp()
    if (.not.check) flag = 58

    check = test_quaternion_log()
    if (.not.check) flag = 59

    check = test_quaternion_pwr()
    if (.not.check) flag = 60

    check = test_quaternion_dot_product()
    if (.not.check) flag = 61

    check = test_quaternion_roll_pitch_yaw()
    if (.not.check) flag = 62

    check = test_line_from_2_points()
    if (.not.check) flag = 63

    check = test_line_from_2_planes()
    if (.not.check) flag = 64

    check = test_line_eval()
    if (.not.check) flag = 65

    check = test_plane_from_3_points()
    if (.not.check) flag = 66

    check = test_plane_from_point_and_normal()
    if (.not.check) flag = 67

    check = test_plane_normal()
    if (.not.check) flag = 68

    check = test_is_parallel_vectors()
    if (.not.check) flag = 69

    check = test_is_parallel_lines()
    if (.not.check) flag = 70

    check = test_is_parallel_planes()
    if (.not.check) flag = 71

    check = test_is_point_on_plane()
    if (.not.check) flag = 72

    check = test_point_to_line_distance()
    if (.not.check) flag = 73

    check = test_is_point_on_line()
    if (.not.check) flag = 74

    check = test_point_to_plane_distance()
    if (.not.check) flag = 75

    check = test_vector_plane_projection()
    if (.not.check) flag = 76

    check = test_fit_line_to_many_points()
    if (.not.check) flag = 77

    check = test_fit_plane_to_many_points()
    if (.not.check) flag = 78

    check = test_flip_plane_normal()
    if (.not.check) flag = 79

    check = test_plucker_line_from_2_points()
    if (.not.check) flag = 80

    check = test_plucker_line_from_line()
    if (.not.check) flag = 81

    check = test_plucker_line_from_2_planes()
    if (.not.check) flag = 82

    check = test_plucker_line_matmul()
    if (.not.check) flag = 83

    check = test_line_from_point_and_vector()
    if (.not.check) flag = 84

    check = test_line_common_normal_1()
    if (.not.check) flag = 85

    check = test_line_common_normal_2()
    if (.not.check) flag = 86

    check = test_line_common_normal_3()
    if (.not.check) flag = 87

    check = test_define_link_csys_1()
    if (.not.check) flag = 88

    check = test_define_link_csys_2()
    if (.not.check) flag = 89

    check = test_define_link_csys_3()
    if (.not.check) flag = 90

    check = test_define_link_csys_4()
    if (.not.check) flag = 91

    check = test_dh_table_1()
    if (.not.check) flag = 92

    check = test_dh_table_2()
    if (.not.check) flag = 93

    check = test_dh_table_3()
    if (.not.check) flag = 94

    check = test_serial_linkage_forward_kinematics()
    if (.not.check) flag = 95

    check = test_serial_linkage_jacobian()
    if (.not.check) flag = 96

    check = test_serial_linkage_inverse_kinematics()
    if (.not.check) flag = 97

    check = test_state_space_initialize()
    if (.not.check) flag = 98

    ! End
    if (flag /= 0) stop flag
end program