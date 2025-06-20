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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module