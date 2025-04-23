module dynamics
    use dynamics_frequency_response
    use dynamics_rotation
    use dynamics_structural
    use dynamics_kinematics
    use dynamics_vibrations
    use dynamics_helper
    use dynamics_stability
    use dynamics_controls
    implicit none
    private

    ! DYNAMICS_FREQUENCY_RESPONSE
    public :: ode_excite
    public :: modal_excite
    public :: harmonic_ode
    public :: frf
    public :: mimo_frf
    public :: chirp
    public :: frequency_response
    public :: frequency_sweep
    public :: compute_modal_damping
    public :: modal_response
    public :: normalize_mode_shapes
    public :: evaluate_accelerance_frf_model
    public :: evaluate_receptance_frf_model
    public :: fit_frf
    public :: FRF_ACCELERANCE_MODEL
    public :: FRF_RECEPTANCE_MODEL
    public :: regression_statistics
    public :: iteration_controls
    public :: lm_solver_options
    public :: convergence_info

    ! DYNAMICS_ROTATION
    public :: rotate_x
    public :: rotate_y
    public :: rotate_z
    public :: rotate
    public :: acceleration_transform
    public :: velocity_transform

    ! DYNAMICS_STRUCTURAL
    public :: DYN_ONE_POINT_INTEGRATION_RULE
    public :: DYN_TWO_POINT_INTEGRATION_RULE
    public :: DYN_THREE_POINT_INTEGRATION_RULE
    public :: DYN_FOUR_POINT_INTEGRATION_RULE
    public :: node
    public :: material
    public :: element
    public :: line_element
    public :: beam_element_2d
    public :: shape_function_derivative
    public :: shape_function_second_derivative
    public :: create_connectivity_matrix
    public :: apply_boundary_conditions
    public :: apply_displacement_constraint
    public :: restore_constrained_values
    public :: point
    public :: beam_element_3d

    ! DYNAMICS_KINEMATICS
    public :: identity_4
    public :: dh_rotate_x
    public :: dh_rotate_z
    public :: dh_translate_x
    public :: dh_translate_z
    public :: dh_matrix
    public :: dh_forward_kinematics
    public :: solve_inverse_kinematics
    public :: vecfcn
    public :: least_squares_solver
    public :: iteration_behavior
    public :: jacobian_generating_vector
    public :: dh_jacobian
    public :: REVOLUTE_JOINT
    public :: PRISMATIC_JOINT

    ! DYNAMICS_VIBRATIONS
    public :: q_factor
    public :: estimate_bandwidth
    public :: logarithmic_decrement
    public :: damping_from_log_decrement
    public :: find_free_response_properties
    public :: rise_time
    public :: find_settling_amplitude
    public :: damping_from_fractional_overshoot
    public :: evaluate_step_response

    ! DYNAMICS_HELPER
    public :: cross_product
    public :: to_skew_symmetric

    ! DYNAMICS_STABILITY
    public :: HYPERBOLIC_FIXED_POINT_SINK
    public :: HYPERBOLIC_FIXED_POINT_SOURCE
    public :: HYPERBOLIC_FIXED_POINT_SADDLE
    public :: NONHYPERBOLIC_FIXED_POINT_UNSTABLE
    public :: NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE
    public :: NONHYPERBOLIC_FIXED_POINT_CENTER
    public :: determine_local_stability

    ! DYNAMICS_CONTROLS
    public :: polynomial
    public :: state_space
    public :: transfer_function
    public :: operator(*)
end module