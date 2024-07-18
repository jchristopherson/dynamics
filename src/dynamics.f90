module dynamics
    use dynamics_frequency_response
    use dynamics_rotation
    implicit none
    private
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
    public :: rotate_x
    public :: rotate_y
    public :: rotate_z

end module