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

    ! End
    if (flag /= 0) stop flag
end program