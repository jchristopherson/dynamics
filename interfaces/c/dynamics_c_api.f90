module dynamics_c_api
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    implicit none

contains
! ******************************************************************************
! DYNAMICS_VIBRATIONS
! ------------------------------------------------------------------------------
function c_q_factor(zeta) result(rst) bind(C, name = "c_q_factor")
    real(c_double), intent(in), value :: zeta
    real(c_double) :: rst
    rst = q_factor(zeta)
end function

! ------------------------------------------------------------------------------
function c_estimate_bandwidth(fn, zeta) result(rst) &
    bind(C, name = "c_estimate_bandwidth")
    real(c_double), intent(in), value :: fn
    real(c_double), intent(in), value :: zeta
    real(c_double) :: rst
    rst = estimate_bandwidth(fn, zeta)
end function

! ------------------------------------------------------------------------------
function c_logarithmic_decrement(x1, x2, n) result(rst) &
    bind(C, name = "c_logarithmic_decrement")
    real(c_double), intent(in), value :: x1
    real(c_double), intent(in), value :: x2
    integer(c_int), intent(in), value :: n
    real(c_double) :: rst
    rst = logarithmic_decrement(x1, x2, n)
end function

! ------------------------------------------------------------------------------
function c_damping_from_log_decrement(delta) result(rst) &
    bind(C, name = "c_damping_from_log_decrement")
    real(c_double), intent(in), value :: delta
    real(c_double) :: rst
    rst = damping_from_log_decrement(delta);
end function

! ------------------------------------------------------------------------------
subroutine c_find_free_response_properties(n, t, x, s, np, delta, fn, x1, x2, &
    t1, t2) bind(C, name = "c_find_free_response_properties")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: t(n)
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in), value :: s
    integer(c_int), intent(in), value :: np
    real(c_double), intent(out) :: delta
    real(c_double), intent(out) :: fn
    real(c_double), intent(out) :: x1
    real(c_double), intent(out) :: x2
    real(c_double), intent(out) :: t1
    real(c_double), intent(out) :: t2
    call find_free_response_properties(t, x, delta, fn, x1 = x1, x2 = x2, &
        t1 = t1, t2 = t2, s = s, n = np)
end subroutine

! ------------------------------------------------------------------------------
function c_rise_time(wn, zeta) result(rst) &
    bind(C, name = "c_rise_time")
    real(c_double), intent(in), value :: wn
    real(c_double), intent(in), value :: zeta
    real(c_double) :: rst
    rst = rise_time(wn, zeta)
end function

! ------------------------------------------------------------------------------
function c_find_settling_amplitude(n, x) result(rst) &
    bind(C, name = "c_find_settling_amplitude")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double) :: rst
    rst = find_settling_amplitude(x)
end function

! ------------------------------------------------------------------------------
function c_damping_from_fractional_overshoot(n, x) result(rst) &
    bind(C, name = "c_damping_from_fractional_overshoot")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double) :: rst
    rst = damping_from_fractional_overshoot(x)
end function

! ------------------------------------------------------------------------------
subroutine c_evaluate_step_response(n, wn, zeta, xs, t, x) &
    bind(C, name = "c_evaluate_step_response")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in), value :: wn
    real(c_double), intent(in), value :: zeta
    real(c_double), intent(in), value :: xs
    real(c_double), intent(in) :: t(n)
    real(c_double), intent(out) :: x(n)
    x = evaluate_step_response(wn, zeta, xs, t)
end subroutine

! ******************************************************************************
! DYNAMICS_ROTATION.F90
! ------------------------------------------------------------------------------
subroutine c_rotate_x(angle, r, ldr) bind(C, name = "c_rotate_x")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (ldr < 3) return
    r(1:3,1:3) = rotate_x(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate_y(angle, r, ldr) bind(C, name = "c_rotate_y")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (ldr < 3) return
    r(1:3,1:3) = rotate_y(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate_z(angle, r, ldr) bind(C, name = "c_rotate_z")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (ldr < 3) return
    r(1:3,1:3) = rotate_z(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate(i, j, k, r, ldr) bind(C, name = "c_rotate")
    real(c_double), intent(in) :: i(3)
    real(c_double), intent(in) :: j(3)
    real(c_double), intent(in) :: k(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,3)
    if (ldr < 3) return
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
    if (ldr < 4) return
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
    if (ldr < 4) return
    r(1:4,1:4) = velocity_transform(omega, v, x)
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module