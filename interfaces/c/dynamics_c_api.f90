module dynamics_c_api
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use ferror
    use diffeq
    use spectrum, only : window
    implicit none

    interface
        subroutine c_vecfcn(nvar, neqn, x, f) bind(C, name = "c_vecfcn")
            use iso_c_binding
            integer(c_int), intent(in), value :: nvar
            integer(c_int), intent(in), value :: neqn
            real(c_double), intent(in) :: x(nvar)
            real(c_double), intent(out) :: f(neqn)
        end subroutine

        subroutine c_modal_excite(n, freq, f) bind(C, name = "c_modal_excite")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: freq
            complex(c_double), intent(out) :: f(n)
        end subroutine

        subroutine c_harmonic_ode(n, freq, t, x, dxdt) &
            bind(C, name = "c_harmonic_ode")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: freq
            real(c_double), intent(in), value :: t
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out) :: dxdt(n)
        end subroutine

        pure function c_window_function(n, bin) result(rst) &
            bind(C, name = "c_window_function")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            integer(c_int), intent(in), value :: bin
            real(c_double) :: rst
        end function
    end interface

    type c_vecfcn_container
        procedure(c_vecfcn), pointer, nopass :: fcn
    end type

    type c_modal_excite_container
        procedure(c_modal_excite), pointer, nopass :: fcn
    end type

    type c_harmonic_ode_container
        procedure(c_harmonic_ode), pointer, nopass :: fcn
    end type

    type, bind(C) :: c_iteration_behavior
        logical(c_bool) :: converge_on_chng
        logical(c_bool) :: converge_on_fcn
        logical(c_bool) :: converge_on_zero_diff
        integer(c_int) :: fcn_count
        integer(c_int) :: gradient_count
        integer(c_int) :: iter_count
        integer(c_int) :: jacobian_count
    end type

    type, bind(C) :: c_frequency_sweep_controls
        integer(c_int) :: cycle_count
        integer(c_int) :: transient_cycles
        integer(c_int) :: points_per_cycle
    end type

    type, bind(C) :: c_iteration_controls
        real(c_double) :: change_in_solution_tolerance
        real(c_double) :: gradient_tolerance
        real(c_double) :: iteration_improvement_tolerance
        real(c_double) :: residual_tolerance
        integer(c_int) :: max_function_evaluations
        integer(c_int) :: max_iteration_between_updates
        integer(c_int) :: max_iteration_count
    end type

    type, bind(C) :: c_regression_statistics
        real(c_double) :: confidence_interval
        real(c_double) :: probability
        real(c_double) :: standard_error
        real(c_double) :: t_statistic
    end type

    type, extends(window) :: c_window
        procedure(c_window_function), pointer, nopass :: fcn
    contains
        procedure, public :: evaluate => cw_eval
    end type

    integer(c_int), parameter :: DYN_RUNGE_KUTTA_23 = 10
    integer(c_int), parameter :: DYN_RUNGE_KUTTA_45 = 11
    integer(c_int), parameter :: DYN_RUNGE_KUTTA_853 = 12
    integer(c_int), parameter :: DYN_ROSENBROCK = 13
    integer(c_int), parameter :: DYN_BDF = 14
    integer(c_int), parameter :: DYN_ADAMS = 15

contains
! ------------------------------------------------------------------------------
subroutine c_report_invalid_input(fcn, name, err)
    character(len = *), intent(in) :: fcn
    character(len = *), intent(in) :: name
    class(errors), intent(inout), optional, target :: err

    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    call errmgr%report_error(fcn, "Invalid Input: " // name, -1)
end subroutine

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
    if (ldr < 3) then
        call c_report_invalid_input("c_rotate_x", "ldr")
        return
    end if
    r(1:3,1:3) = rotate_x(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate_y(angle, r, ldr) bind(C, name = "c_rotate_y")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (ldr < 3) then
        call c_report_invalid_input("c_rotate_y", "ldr")
        return
    end if
    r(1:3,1:3) = rotate_y(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate_z(angle, r, ldr) bind(C, name = "c_rotate_z")
    real(c_double), intent(in), value :: angle
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr, 3)
    if (ldr < 3) then
        call c_report_invalid_input("c_rotate_z", "ldr")
        return
    end if
    r(1:3,1:3) = rotate_z(angle)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_rotate(i, j, k, r, ldr) bind(C, name = "c_rotate")
    real(c_double), intent(in) :: i(3)
    real(c_double), intent(in) :: j(3)
    real(c_double), intent(in) :: k(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,3)
    if (ldr < 3) then
        call c_report_invalid_input("c_rotate", "ldr")
        return
    end if
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
    if (ldr < 4) then
        call c_report_invalid_input("c_acceleration_transform", "ldr")
        return
    end if
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
    if (ldr < 4) then
        call c_report_invalid_input("c_velocity_transform", "ldr")
        return
    end if
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
    
    if (lda < n) then
        call c_report_invalid_input("c_determine_local_stability", "lda")
        return
    end if
    flag = determine_local_stability(a(1:n,1:n), ev = ev)
end subroutine

! ******************************************************************************
! DYNAMICS_KINEMATICS.F90
! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics(n, alpha, a, theta, d, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: alpha(n)
    real(c_double), intent(in) :: a(n)
    real(c_double), intent(in) :: theta(n)
    real(c_double), intent(in) :: d(n)
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics", "ldt")
        return
    end if
    T(1:4,1:4) = dh_forward_kinematics(alpha, a, theta, d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_2(T1, ldt1, T2, ldt2, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_2")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_2", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_2", "ldt2")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_2", "ldt")
        return
    end if
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
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_3", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_3", "ldt2")
        return
    end if
    if (ldt3 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_3", "ldt3")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_3", "ldt")
        return
    end if
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
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_4", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_4", "ldt2")
        return
    end if
    if (ldt3 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_4", "ldt3")
        return
    end if
    if (ldt4 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_4", "ldt4")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_4", "ldt")
        return
    end if
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
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_5", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_5", "ldt2")
        return
    end if
    if (ldt3 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_5", "ldt3")
        return
    end if
    if (ldt4 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_5", "ldt4")
        return
    end if
    if (ldt5 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_5", "ldt5")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_5", "ldt")
        return
    end if
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
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt2")
        return
    end if
    if (ldt3 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt3")
        return
    end if
    if (ldt4 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt4")
        return
    end if
    if (ldt5 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt5")
        return
    end if
    if (ldt6 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt6")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_6", "ldt")
        return
    end if
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
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt2")
        return
    end if
    if (ldt3 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt3")
        return
    end if
    if (ldt4 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt4")
        return
    end if
    if (ldt5 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt5")
        return
    end if
    if (ldt6 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt6")
        return
    end if
    if (ldt7 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt7")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_7", "ldt")
        return
    end if
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
    if (ldt1 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt1")
        return
    end if
    if (ldt2 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt2")
        return
    end if
    if (ldt3 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt3")
        return
    end if
    if (ldt4 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt4")
        return
    end if
    if (ldt5 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt5")
        return
    end if
    if (ldt6 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt6")
        return
    end if
    if (ldt7 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt7")
        return
    end if
    if (ldt8 < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt8")
        return
    end if
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_forward_kinematics_8", "ldt")
        return
    end if
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
    if (ldjac < 6) then
        call c_report_invalid_input("c_dh_jacobian", "ldjac")
        return
    end if
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
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_matrix", "ldt")
        return
    end if
    T(1:4,1:4) = dh_matrix(alpha, a, theta, d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_rotate_x(alpha, T, ldt) &
    bind(C, name = "c_dh_rotate_x")
    real(c_double), intent(in), value :: alpha
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_rotate_x", "ldt")
        return
    end if
    T(1:4,1:4) = dh_rotate_x(alpha)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_rotate_z(theta, T, ldt) &
    bind(C, name = "c_dh_rotate_z")
    real(c_double), intent(in), value :: theta
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_rotate_z", "ldt")
        return
    end if
    T(1:4,1:4) = dh_rotate_z(theta)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_translate_x(a, T, ldt) &
    bind(C, name = "c_dh_translate_x")
    real(c_double), intent(in), value :: a
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_translate_x", "ldt")
        return
    end if
    T(1:4,1:4) = dh_translate_x(a)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_translate_z(d, T, ldt) &
    bind(C, name = "c_dh_translate_z")
    real(c_double), intent(in), value :: d
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) then
        call c_report_invalid_input("c_dh_translate_z", "ldt")
        return
    end if
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
    if (ldr < 3) then
        call c_report_invalid_input("c_jacobian_generating_vector", "ldr")
        return
    end if
    jvec = jacobian_generating_vector(d, k, R(1:3,1:3), jtype)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_solve_inverse_kinematics(njoints, neqn, mdl, qo, constraints, &
    jvar, resid, ib) bind(C, name = "c_solve_inverse_kinematics")
    integer(c_int), intent(in), value :: njoints
    integer(c_int), intent(in), value :: neqn
    type(c_funptr), intent(in), value :: mdl
    real(c_double), intent(in) :: qo(njoints)
    real(c_double), intent(in) :: constraints(neqn)
    real(c_double), intent(out) :: jvar(njoints)
    real(c_double), intent(out) :: resid(neqn)
    type(c_iteration_behavior), intent(out) :: ib
    type(iteration_behavior) :: iter
    procedure(vecfcn), pointer :: fcn
    procedure(c_vecfcn), pointer :: fptr
    type(c_vecfcn_container) :: arg
    call c_f_procpointer(mdl, fptr)
    fcn => sik_fcn
    arg%fcn => fptr
    
    jvar = solve_inverse_kinematics(fcn, qo, constraints, df = resid, &
        ib = iter, args = arg)
    ib%converge_on_chng = iter%converge_on_chng
    ib%converge_on_fcn = iter%converge_on_fcn
    ib%converge_on_zero_diff = iter%converge_on_zero_diff
    ib%fcn_count = iter%fcn_count
    ib%gradient_count = iter%gradient_count
    ib%iter_count = iter%iter_count
    ib%jacobian_count = iter%jacobian_count
end subroutine

! --------------------
subroutine sik_fcn(x, f, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args
    select type (args)
    class is (c_vecfcn_container)
        call args%fcn(size(x), size(f), x, f)
    end select
end subroutine

! ******************************************************************************
! DYNAMICS_FREQUENCY_RESPONSE.F90
! ------------------------------------------------------------------------------
subroutine c_frequency_response(n, nfreq, mass, ldm, stiff, ldk, alpha, beta, &
    freq, frc, modes, modeshapes, ldms, rsp, ldr) &
    bind(C, name = "c_frequency_response")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: nfreq
    integer(c_int), intent(in), value :: ldm
    integer(c_int), intent(in), value :: ldk
    integer(c_int), intent(in), value :: ldms
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(in) :: mass(ldm,n)
    real(c_double), intent(in) :: stiff(ldk,n)
    real(c_double), intent(in), value :: alpha
    real(c_double), intent(in), value :: beta
    real(c_double), intent(in) :: freq(nfreq)
    type(c_funptr), intent(in), value :: frc
    real(real64), intent(out) :: modes(n)
    real(real64), intent(out) :: modeshapes(ldms,n)
    complex(real64), intent(out) :: rsp(ldr,n)

    type(frf) :: frsp
    type(c_modal_excite_container) :: arg
    procedure(c_modal_excite), pointer :: fptr
    procedure(modal_excite), pointer :: fcn
    real(real64), allocatable, dimension(:) :: mds
    real(real64), allocatable, dimension(:,:) :: ms

    if (ldm < n) then
        call c_report_invalid_input("c_frequency_response", "ldm")
        return
    end if
    if (ldk < n) then
        call c_report_invalid_input("c_frequency_response", "ldk")
        return
    end if
    if (ldms < n) then
        call c_report_invalid_input("c_frequency_response", "ldms")
        return
    end if
    if (ldr < nfreq) then
        call c_report_invalid_input("c_frequency_response", "ldr")
        return
    end if

    call c_f_procpointer(frc, fptr)
    arg%fcn => fptr
    fcn => cfr_fcn
    
    frsp = frequency_response(mass(1:n,1:n), stiff(1:n,1:n), alpha, beta, &
        freq, fcn, modes = mds, modeshapes = ms, args = arg)
    rsp(1:nfreq,1:n) = frsp%responses
    modes = mds
    modeshapes(1:n,1:n) = ms
end subroutine

! --------------------
subroutine cfr_fcn(freq, frc, args)
    real(real64), intent(in) :: freq
    complex(real64), intent(out), dimension(:) :: frc
    class(*), intent(inout), optional :: args
    select type (args)
    class is (c_modal_excite_container)
        call args%fcn(size(frc), freq, frc)
    end select
end subroutine

! ------------------------------------------------------------------------------
function c_compute_modal_damping(lambda, alpha, beta) result(rst) &
    bind(C, name = "c_compute_modal_damping")
    real(c_double), intent(in), value :: lambda
    real(c_double), intent(in), value :: alpha
    real(c_double), intent(in), value :: beta
    real(c_double) :: rst
    rst = compute_modal_damping(lambda, alpha, beta)
end function

! ------------------------------------------------------------------------------
function c_chirp(t, amp, span, f1Hz, f2Hz) result(rst) bind(C, name = "c_chirp")
    real(c_double), intent(in), value :: t
    real(c_double), intent(in), value :: amp
    real(c_double), intent(in), value :: span
    real(c_double), intent(in), value :: f1Hz
    real(c_double), intent(in), value :: f2Hz
    real(c_double) :: rst
    rst = chirp(t, amp, span, f1Hz, f2Hz)
end function

! ------------------------------------------------------------------------------
subroutine c_modal_response(n, mass, ldm, stiff, ldk, freqs, modeshapes, ldms) &
    bind(C, name = "c_modal_response")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: ldm
    integer(c_int), intent(in), value :: ldk
    integer(c_int), intent(in), value :: ldms
    real(c_double), intent(in) :: mass(ldm,n)
    real(c_double), intent(in) :: stiff(ldk,n)
    real(c_double), intent(out) :: freqs(n)
    real(c_double), intent(out) :: modeshapes(ldms,n)

    real(real64), allocatable, dimension(:) :: mds
    real(real64), allocatable, dimension(:,:) :: ms

    if (ldm < n) then
        call c_report_invalid_input("c_modal_response", "ldm")
        return
    end if
    if (ldk < n) then
        call c_report_invalid_input("c_modal_response", "ldk")
        return
    end if
    if (ldms < n) then
        call c_report_invalid_input("c_modal_response", "ldms")
        return
    end if

    
    call modal_response(mass(1:n,1:n), stiff(1:n,1:n), mds, ms)
    freqs = mds
    modeshapes(1:n,1:n) = ms
end subroutine

! ------------------------------------------------------------------------------
subroutine c_normalize_mode_shapes(n, x, ldx) &
    bind(C, name = "c_normalize_mode_shapes")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: ldx
    real(c_double), intent(inout) :: x(ldx,n)
    if (ldx < n) then
        call c_report_invalid_input("c_normalize_mode_shapes", "ldx")
        return
    end if
    call normalize_mode_shapes(x(1:n,1:n))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_frf_sweep(n, nfreq, fcn, freq, iv, solver, rsp, ldr, opts) &
    bind(C, name = "c_frf_sweep")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: nfreq
    integer(c_int), intent(in), value :: ldr
    type(c_funptr), intent(in), value :: fcn
    real(c_double), intent(in) :: freq(nfreq)
    real(c_double), intent(in) :: iv(n)
    integer(c_int), intent(in), value :: solver
    complex(c_double), intent(out) :: rsp(ldr,n)
    type(c_frequency_sweep_controls), intent(in) :: opts

    type(c_harmonic_ode_container) :: arg
    procedure(c_harmonic_ode), pointer :: fptr
    procedure(harmonic_ode), pointer :: odefcn
    type(runge_kutta_23), target :: rk23
    type(runge_kutta_45), target :: rk45
    type(runge_kutta_853), target :: rk853
    type(rosenbrock), target :: rbrk
    type(bdf), target :: bdiff
    type(adams), target :: pece
    class(ode_integrator), pointer :: integrator

    type(frf) :: frsp

    if (ldr < nfreq) then
        call c_report_invalid_input("c_frf_sweep", "ldr")
        return
    end if

    call c_f_procpointer(fcn, fptr)
    arg%fcn => fptr
    odefcn => cfrf_sweep_fcn
    

    select case (solver)
    case (DYN_ADAMS)
        integrator => pece
    case (DYN_BDF)
        integrator => bdiff
    case (DYN_ROSENBROCK)
        integrator => rbrk
    case (DYN_RUNGE_KUTTA_23)
        integrator => rk23
    case (DYN_RUNGE_KUTTA_45)
        integrator => rk45
    case (DYN_RUNGE_KUTTA_853)
        integrator => rk853
    case default
        integrator => rk45
    end select

    frsp = frequency_sweep(odefcn, freq, iv, solver = integrator, args = arg, &
        ncycles = opts%cycle_count, ntransient = opts%transient_cycles, &
        points = opts%points_per_cycle)
    rsp(1:nfreq,1:n) = frsp%responses
end subroutine

! --------------------
subroutine cfrf_sweep_fcn(freq, t, x, dxdt, args)
    real(real64), intent(in) :: freq, t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: dxdt
    class(*), intent(inout), optional :: args

    select type (args)
    class is (c_harmonic_ode_container)
        call args%fcn(size(x), freq, t, x, dxdt)
    end select
end subroutine

! ------------------------------------------------------------------------------
subroutine c_set_frequency_sweep_defaults(x) &
    bind(C, name = "c_set_frequency_sweep_defaults")
    type(c_frequency_sweep_controls), intent(inout) :: x
    x%cycle_count = 20
    x%transient_cycles = 200
    x%points_per_cycle = 1000
end subroutine

! ------------------------------------------------------------------------------
subroutine c_evaluate_accelerance_frf_model(n, norder, mdl, omega, h) &
    bind(C, name = "c_evaluate_accelerance_frf_model")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: norder
    real(c_double), intent(in) :: mdl(3 * norder)
    real(c_double), intent(in) :: omega(n)
    complex(c_double), intent(out) :: h(n)

    h = evaluate_accelerance_frf_model(mdl, omega)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_evaluate_receptance_frf_model(n, norder, mdl, omega, h) &
    bind(C, name = "c_evaluate_receptance_frf_model")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: norder
    real(c_double), intent(in) :: mdl(3 * norder)
    real(c_double), intent(in) :: omega(n)
    complex(c_double), intent(out) :: h(n)

    h = evaluate_receptance_frf_model(mdl, omega)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_set_iteration_controls_defaults(x) &
    bind(C, name = "c_set_iteration_controls_defaults")
    type(c_iteration_controls), intent(inout) :: x
    type(iteration_controls) :: c
    call c%set_to_default()
    x%change_in_solution_tolerance = c%change_in_solution_tolerance
    x%gradient_tolerance = c%gradient_tolerance
    x%iteration_improvement_tolerance = c%iteration_improvement_tolerance
    x%max_function_evaluations = c%max_function_evaluations
    x%max_iteration_between_updates = c%max_iteration_between_updates
    x%max_iteration_count = c%max_iteration_count
    x%residual_tolerance = c%residual_tolerance
end subroutine

! ------------------------------------------------------------------------------
subroutine c_fit_frf(n, norder, method, freq, rsp, maxp, minp, controls, mdl, &
    stats) bind(C, name = "c_fit_frf")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: norder
    integer(c_int), intent(in), value :: method
    real(c_double), intent(in) :: freq(n)
    complex(c_double), intent(in) :: rsp(n)
    real(c_double), intent(in) :: maxp(3 * norder)
    real(c_double), intent(in) :: minp(3 * norder)
    type(c_iteration_controls), intent(in) :: controls
    real(c_double), intent(out) :: mdl(3 * norder)
    type(c_regression_statistics), intent(out) :: stats(3 * norder)

    integer(int32) :: i
    type(iteration_controls) :: cntrls
    type(regression_statistics) :: fs(3 * norder)

    

    cntrls%change_in_solution_tolerance = controls%change_in_solution_tolerance
    cntrls%gradient_tolerance = controls%gradient_tolerance
    cntrls%iteration_improvement_tolerance = controls%iteration_improvement_tolerance
    cntrls%max_function_evaluations = controls%max_function_evaluations
    cntrls%max_iteration_between_updates = controls%max_iteration_between_updates
    cntrls%max_iteration_count = controls%max_iteration_count
    cntrls%residual_tolerance = controls%residual_tolerance

    mdl = fit_frf(method, norder, freq, rsp, maxp = maxp, minp = minp, &
        stats = fs, controls = cntrls)
    do i = 1, size(fs)
        stats(i)%confidence_interval = fs(i)%confidence_interval
        stats(i)%probability = fs(i)%probability
        stats(i)%standard_error = fs(i)%standard_error
        stats(i)%t_statistic = fs(i)%t_statistic
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function cw_eval(this, bin) result(rst)
    class(c_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    rst = this%fcn(this%size, bin)
end function

! ------------------------------------------------------------------------------
subroutine c_siso_frequency_response(n, nf, x, y, fs, winsize, winfun, method, &
    freq, rsp) bind(C, name = "c_siso_frequency_response")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: nf
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(in), value :: fs
    integer(c_int), intent(in), value :: winsize
    type(c_funptr), intent(in), value :: winfun
    integer(c_int), intent(in), value :: method
    real(c_double), intent(out) :: freq(nf)
    complex(c_double), intent(out) :: rsp(nf)

    type(c_window) :: win
    type(frf) :: frsp
    procedure(c_window_function), pointer :: fcn
    integer(c_int) :: m

    
    if (mod(winsize, 2) == 0) then
        m = winsize / 2 + 1
    else
        m = (winsize + 1) / 2
    end if
    if (nf /= m) return

    call c_f_procpointer(winfun, fcn)
    win%size = winsize
    win%fcn => fcn

    frsp = frequency_response(x, y, fs, win = win, method = method)
    freq = frsp%frequency
    rsp = frsp%responses(:,1)
end subroutine

! ******************************************************************************
! DYNAMICS_HELPER.F90
! ------------------------------------------------------------------------------
subroutine c_cross_product(x, y, z) bind(C, name = "c_cross_product")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double), intent(out) :: z(3)
    z = cross_product(x, y)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_to_skew_symmetric(x, y, ldy) bind(C, name = "c_to_skew_symmetric")
    real(c_double), intent(in) :: x(3)
    integer(c_int), intent(in), value :: ldy
    real(c_double), intent(out) :: y(ldy,3)
    if (ldy < 3) then
        call c_report_invalid_input("c_to_skew_symmetric", "ldy")
        return
    end if
    y(1:3,1:3) = to_skew_symmetric(x)
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module