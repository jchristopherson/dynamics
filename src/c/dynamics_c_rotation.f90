module dynamics_c_rotation
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use diffeq
    use spectrum, only : window
    use dynamics_error_handling
    use nonlin
    use dynamics_c_types
    implicit none

contains

subroutine c_rotate_x(angle, r, ldr) bind(C, name = "c_rotate_x")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (c_api_error(ldr < 3, DYN_INVALID_INPUT_ERROR, &
        "c_rotate_x: ldr must be >= 3.")) return
    r(1:3,1:3) = rotate_x(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate_y(angle, r, ldr) bind(C, name = "c_rotate_y")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (c_api_error(ldr < 3, DYN_INVALID_INPUT_ERROR, &
        "c_rotate_y: ldr must be >= 3.")) return
    r(1:3,1:3) = rotate_y(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate_z(angle, r, ldr) bind(C, name = "c_rotate_z")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (c_api_error(ldr < 3, DYN_INVALID_INPUT_ERROR, &
        "c_rotate_z: ldr must be >= 3.")) return
    r(1:3,1:3) = rotate_z(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate(i, j, k, r, ldr) bind(C, name = "c_rotate")
    real(c_double), intent(in) :: i(3)
    real(c_double), intent(in) :: j(3)
    real(c_double), intent(in) :: k(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,3)
    if (c_api_error(ldr < 3, DYN_INVALID_INPUT_ERROR, &
        "c_rotate: ldr must be >= 3.")) return
    r(1:3,1:3) = rotate(i, j, k)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_acceleration_transform(alpha, omega, a, x, r, ldr) &
    bind(C, name = "c_acceleration_transform")
    real(c_double), intent(in) :: alpha(3)
    real(c_double), intent(in) :: omega(3)
    real(c_double), intent(in) :: a(3)
    real(c_double), intent(in) :: x(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,4)
    if (c_api_error(ldr < 4, DYN_INVALID_INPUT_ERROR, &
        "c_acceleration_transform: ldr must be >= 4.")) return
    r(1:4,1:4) = acceleration_transform(alpha, omega, a, x)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_velocity_transform(omega, v, x, r, ldr) &
    bind(C, name = "c_velocity_transform")
    real(c_double), intent(in) :: omega(3)
    real(c_double), intent(in) :: v(3)
    real(c_double), intent(in) :: x(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,4)
    if (c_api_error(ldr < 4, DYN_INVALID_INPUT_ERROR, &
        "c_velocity_transform: ldr must be >= 4.")) return
    r(1:4,1:4) = velocity_transform(omega, v, x)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_to_angle_axis(r, ldr, angle, axis) bind(C, name = "c_to_angle_axis")
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(in) :: r(ldr,3)
    real(c_double), intent(out) :: angle
    real(c_double), intent(out) :: axis(3)
    call to_angle_axis(r(1:3,1:3), angle, axis)
end subroutine

end module
