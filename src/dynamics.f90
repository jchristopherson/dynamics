module dynamics
    use dynamics_frequency_response
    use dynamics_rotation
    use dynamics_structural
    use dynamics_kinematics
    implicit none
    private
    ! DYNAMICS_FREQUENCY_RESPONSE
    public :: ode_excite
    public :: modal_excite
    public :: harmonic_ode
    public :: frf
    public :: chirp
    public :: frequency_response
    public :: frequency_sweep
    public :: compute_modal_damping
    public :: modal_response
    public :: normalize_mode_shapes

    ! DYNAMICS_ROTATION
    public :: rotate_x
    public :: rotate_y
    public :: rotate_z

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

    ! DYNAMICS_KINEMATICS
    public :: identity_4
    public :: dh_rotate_x
    public :: dh_rotate_z
    public :: dh_translate_x
    public :: dh_translate_z
    public :: dh_matrix
    public :: dh_forward_kinematics
end module