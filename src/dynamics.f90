module dynamics
    use dynamics_frequency_response
    use dynamics_rotation
    use dynamics_structural
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
end module