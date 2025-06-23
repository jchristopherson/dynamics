module dynamics_c_api
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use ferror
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

! ******************************************************************************
! DYNAMICS_STABILITY.F90
! ------------------------------------------------------------------------------
subroutine c_determine_local_stability(n, a, lda, ev, flag) &
    bind(C, name = "c_determine_local_stability")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(in) :: a(lda,n)
    complex(c_double), intent(out) :: ev(n)
    integer(c_int), intent(out) :: flag
    type(errors) :: err
    if (lda < n) return;
    call err%set_exit_on_error(.false.)
    flag = determine_local_stability(a(1:n,1:n), ev = ev, err = err)
end subroutine

! ******************************************************************************
! DYNAMICS_KINEMATICS.F90
! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_2(T1, ldt1, T2, ldt2, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_2")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_3(T1, ldt1, T2, ldt2, T3, ldt3, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_3")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt3 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_4(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_4")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt3 < 4) return;
    if (ldt4 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_5(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_5")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt3 < 4) return;
    if (ldt4 < 4) return;
    if (ldt5 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_6(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T6, ldt6, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_6")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt6
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(in) :: T6(ldt6,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt3 < 4) return;
    if (ldt4 < 4) return;
    if (ldt5 < 4) return;
    if (ldt6 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4), T6(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_7(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T6, ldt6, T7, ldt7, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_7")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt6
    integer(c_int), intent(in), value :: ldt7
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(in) :: T6(ldt6,4)
    real(c_double), intent(in) :: T7(ldt7,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt3 < 4) return;
    if (ldt4 < 4) return;
    if (ldt5 < 4) return;
    if (ldt6 < 4) return;
    if (ldt7 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4), T6(1:4,1:4), T7(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_8(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T6, ldt6, T7, ldt7, T8, ldt8, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_8")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt6
    integer(c_int), intent(in), value :: ldt7
    integer(c_int), intent(in), value :: ldt8
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(in) :: T6(ldt6,4)
    real(c_double), intent(in) :: T7(ldt7,4)
    real(c_double), intent(in) :: T8(ldt8,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) return;
    if (ldt2 < 4) return;
    if (ldt3 < 4) return;
    if (ldt4 < 4) return;
    if (ldt5 < 4) return;
    if (ldt6 < 4) return;
    if (ldt7 < 4) return;
    if (ldt8 < 4) return;
    if (ldt < 4) return;
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4), T6(1:4,1:4), T7(1:4,1:4), T8(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_jacobian(n, alpha, a, theta, d, jtypes, jac, ldjac) &
    bind(C, name = "c_dh_jacobian")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: ldjac
    real(c_double), intent(in) :: alpha(n)
    real(c_double), intent(in) :: a(n)
    real(c_double), intent(in) :: theta(n)
    real(c_double), intent(in) :: d(n)
    integer(c_int), intent(in) :: jtypes(n)
    real(c_double), intent(out) :: jac(ldjac,n)
    if (ldjac < 6) return
    jac(1:6,1:n) = dh_jacobian(alpha, a, theta, d, jtypes)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_matrix(alpha, a, theta, d, T, ldt) &
    bind(C, name = "c_dh_matrix")
    real(c_double), intent(in), value :: alpha
    real(c_double), intent(in), value :: a
    real(c_double), intent(in), value :: theta
    real(c_double), intent(in), value :: d
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) return
    T(1:4,1:4) = dh_matrix(alpha, a, theta, d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_rotate_x(alpha, T, ldt) &
    bind(C, name = "c_dh_rotate_x")
    real(c_double), intent(in), value :: alpha
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) return
    T(1:4,1:4) = dh_rotate_x(alpha)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_rotate_z(theta, T, ldt) &
    bind(C, name = "c_dh_rotate_z")
    real(c_double), intent(in), value :: theta
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) return
    T(1:4,1:4) = dh_rotate_z(theta)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_translate_x(a, T, ldt) &
    bind(C, name = "c_dh_translate_x")
    real(c_double), intent(in), value :: a
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) return
    T(1:4,1:4) = dh_translate_x(a)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_translate_z(d, T, ldt) &
    bind(C, name = "c_dh_translate_z")
    real(c_double), intent(in), value :: d
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) return
    T(1:4,1:4) = dh_translate_z(d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_jacobian_generating_vector(d, k, R, ldr, jtype, jvec) &
    bind(C, name = "c_jacobian_generating_vector")
    real(c_double), intent(in) :: d(3)
    real(c_double), intent(in) :: k(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(in) :: R(ldr,3)
    integer(c_int), intent(in), value :: jtype
    real(c_double), intent(out) :: jvec(6)
    if (ldr < 3) return
    jvec = jacobian_generating_vector(d, k, R(1:3,1:3), jtype)
end subroutine

! ------------------------------------------------------------------------------
end module