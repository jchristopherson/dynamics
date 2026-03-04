module dynamics_c_api
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use ferror
    use diffeq
    use spectrum, only : window
    use dynamics_quaternions
    use dynamics_geometry
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

        subroutine c_constraint_equations(n, neqn, nparam, xg, fg, xc, p, &
            fc) bind(C, name = "c_constraint_equations")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            integer(c_int), intent(in), value :: neqn
            integer(c_int), intent(in), value :: nparam
            real(c_double), intent(in) :: xg(n)
            real(c_double), intent(in) :: fg(n)
            real(c_double), intent(in) :: xc(neqn)
            real(c_double), intent(in) :: p(nparam)
            real(c_double), intent(out) :: fc(neqn)
        end subroutine

        subroutine c_ode_fit(neqn, nparam, mdl, t, x, F, dxdt) &
            bind(C, name = "c_ode_fit")
            use iso_c_binding
            integer(c_int), intent(in), value :: neqn
            integer(c_int), intent(in), value :: nparam
            real(c_double), intent(in) :: mdl(nparam)
            real(c_double), intent(in), value :: t
            real(c_double), intent(in) :: x(neqn)
            real(c_double), intent(in), value :: F
            real(c_double), intent(out) :: dxdt(neqn)
        end subroutine
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

    type c_siso_fit_container
        procedure(c_ode_fit), pointer, nopass :: odefcn
        procedure(c_constraint_equations), pointer, nopass :: constraints
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
    
    type, bind(C) :: c_lm_solver_options
        real(c_double) :: damping_decrease_factor
        real(c_double) :: damping_increase_factor
        real(c_double) :: finite_difference_step_size
        integer(c_int) :: method
    end type

    type, extends(window) :: c_window
        procedure(c_window_function), pointer, nopass :: fcn
    contains
        procedure, public :: evaluate => cw_eval
    end type

    type, bind(C) :: c_dynamic_system_measurement
        integer(c_int) :: npts
        type(c_ptr) :: input
        type(c_ptr) :: output
        type(c_ptr) :: t
    end type

    integer(c_int), parameter :: DYN_RUNGE_KUTTA_23 = 10
    integer(c_int), parameter :: DYN_RUNGE_KUTTA_45 = 11
    integer(c_int), parameter :: DYN_RUNGE_KUTTA_853 = 12
    integer(c_int), parameter :: DYN_ROSENBROCK = 13
    integer(c_int), parameter :: DYN_BDF = 14
    integer(c_int), parameter :: DYN_ADAMS = 15

    type, bind(C) :: c_quaternion
        real(c_double) :: w
        real(c_double) :: x
        real(c_double) :: y
        real(c_double) :: z
    end type

    type, bind(C) :: c_plane
        real(c_double) :: a
        real(c_double) :: b
        real(c_double) :: c
        real(c_double) :: d
    end type

    type, bind(C) :: c_line
        real(c_double) :: r0(3)
        real(c_double) :: v(3)
    end type

    type, bind(C) :: c_plucker_line
        real(c_double) :: u(3);
        real(c_double) :: m(3);
    end type

    interface assignment(=)
        module procedure :: convert_to_c_quaternion
        module procedure :: convert_from_c_quaternion
        module procedure :: convert_to_c_line
        module procedure :: convert_from_c_line
        module procedure :: convert_to_c_plane
        module procedure :: convert_from_c_plane
        module procedure :: convert_to_c_plucker_line
        module procedure :: convert_from_c_plucker_line
    end interface

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

! ------------------------------------------------------------------------------
subroutine c_to_angle_axis(r, ldr, angle, axis) bind(C, name = "c_to_angle_axis")
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(in) :: r(ldr,3)
    real(c_double), intent(out) :: angle
    real(c_double), intent(out) :: axis(3)
    call to_angle_axis(r(1:3,1:3), angle, axis)
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
    class(ode_integrator), pointer :: integrator_obj

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
        integrator_obj => pece
    case (DYN_BDF)
        integrator_obj => bdiff
    case (DYN_ROSENBROCK)
        integrator_obj => rbrk
    case (DYN_RUNGE_KUTTA_23)
        integrator_obj => rk23
    case (DYN_RUNGE_KUTTA_45)
        integrator_obj => rk45
    case (DYN_RUNGE_KUTTA_853)
        integrator_obj => rk853
    case default
        integrator_obj => rk45
    end select

    frsp = frequency_sweep(odefcn, freq, iv, solver = integrator_obj, &
        args = arg, ncycles = opts%cycle_count, &
        ntransient = opts%transient_cycles, points = opts%points_per_cycle)
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
function c_vector_angle(x, y) result(rst) bind(C, name = "c_vector_angle")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double) :: rst
    rst = vector_angle(x, y)
end function

! ------------------------------------------------------------------------------
function c_scalar_projection(x, y) result(rst) bind(C, name = "c_scalar_projection")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double) :: rst
    rst = scalar_projection(x, y)
end function

! ------------------------------------------------------------------------------
subroutine c_vector_projection(x, y, z) bind(C, name = "c_vector_projection")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double), intent(out) :: z(3)
    z = vector_projection(x, y)
end subroutine

! ------------------------------------------------------------------------------
function c_vector_magnitude(n, x) result(rst) bind(C, name = "c_vector_magnitude")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double) :: rst
    rst = norm2(x)
end function

! ------------------------------------------------------------------------------
subroutine c_vector_normalize(n, x) bind(C, name = "c_vector_normalize")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: x(n)
    x = x / norm2(x)
end subroutine

! ------------------------------------------------------------------------------
function c_dot_product(n, x, y) result(rst) bind(C, name = "c_dot_product")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double) :: rst
    rst = dot_product(x, y)
end function

! ******************************************************************************
! DYNAMICS_SYSTEM_ID.F90
! ------------------------------------------------------------------------------
subroutine c_siso_model_fit_least_squares(nsets, nparams, neqns, fcn, x, ic, &
    p, integrator, ind, maxp, minp, controls, opts, nconstraints, xc, yc, &
    constraints, nweights, weights, stats, info) &
    bind(C, name = "c_siso_model_fit_least_squares")
    integer(c_int), intent(in), value :: nsets
    integer(c_int), intent(in), value :: nparams
    integer(c_int), intent(in), value :: neqns
    type(c_funptr), intent(in), value :: fcn
    type(c_dynamic_system_measurement), intent(in) :: x(nsets)
    real(c_double), intent(in) :: ic(neqns)
    real(c_double), intent(inout) :: p(nparams)
    integer(c_int), intent(in), value :: integrator
    integer(c_int), intent(in), value :: ind
    real(c_double), intent(in) :: maxp(nparams)
    real(c_double), intent(in) :: minp(nparams)
    type(c_iteration_controls), intent(in) :: controls
    type(c_lm_solver_options), intent(in) :: opts
    integer(c_int), intent(in), value :: nconstraints
    real(c_double), intent(in) :: xc(nconstraints)
    real(c_double), intent(in) :: yc(nconstraints)
    type(c_funptr), intent(in), value :: constraints
    integer(c_int), intent(in), value :: nweights
    real(c_double), intent(in) :: weights(nweights)
    type(c_regression_statistics), intent(out) :: stats(nparams)
    type(c_iteration_behavior), intent(out) :: info

    ! Variables
    logical :: uses_constraints, uses_weights
    integer(int32) :: i, nw, flag(1)
    type(dynamic_system_measurement), allocatable, dimension(:) :: fx
    procedure(c_ode_fit), pointer :: f_ode
    procedure(c_constraint_equations), pointer :: f_constraints
    procedure(ode), pointer :: odeptr
    procedure(constraint_equations), pointer :: constraints_pointer
    type(c_siso_fit_container) :: args
    real(real64), pointer, dimension(:) :: temp
    type(regression_statistics), allocatable, dimension(:) :: f_stats
    type(convergence_info) :: f_info
    type(iteration_controls) :: f_controls
    type(runge_kutta_23), target :: rk23
    type(runge_kutta_45), target :: rk45
    type(runge_kutta_853), target :: rk853
    type(rosenbrock), target :: rbrk
    type(bdf), target :: bdiff
    type(adams), target :: pece
    class(ode_integrator), pointer :: integrator_obj
    type(lm_solver_options) :: f_opt

    ! Uses constraints?
    if (nconstraints == 0 .or. .not.c_associated(constraints)) then
        uses_constraints = .false.
    else
        uses_constraints = .true.
    end if

    ! Establish function pointers
    call c_f_procpointer(fcn, f_ode)
    args%odefcn => f_ode
    odeptr => siso_fit_ode
    if (uses_constraints) then
        call c_f_procpointer(constraints, f_constraints)
        args%constraints => f_constraints
        constraints_pointer => siso_constraint_equations
    end if

    ! Uses weights?
    if (nweights == 0) then
        uses_weights = .false.
    else
        uses_weights = .true.
        nw = 0
        do i = 1, nsets
            nw = nw + x(i)%npts
        end do
        if (nweights /= nw) then
            call c_report_invalid_input("c_siso_model_fit_least_squares", &
                "nweights")
            return
        end if
    end if

    ! Convert the inputs
    allocate(fx(nsets))
    do i = 1, nsets
        flag(1) = x(i)%npts
        call c_f_pointer(x(i)%input, temp, flag)
        allocate(fx(i)%input(x(i)%npts), source = temp)
        
        call c_f_pointer(x(i)%output, temp, flag)
        allocate(fx(i)%output(x(i)%npts), source = temp)

        call c_f_pointer(x(i)%t, temp, flag)
        allocate(fx(i)%t(x(i)%npts), source = temp)
    end do

    ! Define the integrator
    select case (integrator)
    case (DYN_ADAMS)
        integrator_obj => pece
    case (DYN_BDF)
        integrator_obj => bdiff
    case (DYN_ROSENBROCK)
        integrator_obj => rbrk
    case (DYN_RUNGE_KUTTA_23)
        integrator_obj => rk23
    case (DYN_RUNGE_KUTTA_45)
        integrator_obj => rk45
    case (DYN_RUNGE_KUTTA_853)
        integrator_obj => rk853
    case default
        integrator_obj => rk45
    end select

    ! Set up the iteration controls
    f_controls%change_in_solution_tolerance = &
        controls%change_in_solution_tolerance
    f_controls%gradient_tolerance = controls%gradient_tolerance
    f_controls%iteration_improvement_tolerance = &
        controls%iteration_improvement_tolerance
    f_controls%max_function_evaluations = controls%max_function_evaluations
    f_controls%max_iteration_between_updates = &
        controls%max_iteration_between_updates
    f_controls%max_iteration_count = controls%max_iteration_count
    f_controls%residual_tolerance = controls%residual_tolerance

    ! Set up the solver options
    f_opt%damping_decrease_factor = opts%damping_decrease_factor
    f_opt%damping_increase_factor = opts%damping_increase_factor
    f_opt%finite_difference_step_size = opts%finite_difference_step_size
    f_opt%method = opts%method

    ! Process
    allocate(f_stats(nparams))
    if (uses_constraints .and. uses_weights) then
        call siso_model_fit_least_squares(odeptr, fx, ic, p, &
            integrator = integrator_obj, ind = ind, maxp = maxp, minp = minp, &
            stats = f_stats, controls = f_controls, info = f_info, xc = xc, &
            yc = yc, constraints = constraints_pointer, weights = weights, &
            args = args, settings = f_opt)
    else if (uses_constraints .and. .not. uses_weights) then
        call siso_model_fit_least_squares(odeptr, fx, ic, p, &
            integrator = integrator_obj, ind = ind, maxp = maxp, minp = minp, &
            stats = f_stats, controls = f_controls, info = f_info, xc = xc, &
            yc = yc, constraints = constraints_pointer, args = args, &
            settings = f_opt)
    else if (.not. uses_constraints .and. uses_weights) then
        call siso_model_fit_least_squares(odeptr, fx, ic, p, &
            integrator = integrator_obj, ind = ind, maxp = maxp, minp = minp, &
            stats = f_stats, controls = f_controls, info = f_info, &
            weights = weights, args = args, settings = f_opt)
    else
        call siso_model_fit_least_squares(odeptr, fx, ic, p, &
            integrator = integrator_obj, ind = ind, maxp = maxp, minp = minp, &
            stats = f_stats, controls = f_controls, info = f_info, &
            args = args, settings = f_opt)
    end if

    ! Extract the output information
    info%converge_on_chng = logical(f_info%converge_on_solution_change, c_bool)
    info%converge_on_fcn = logical(f_info%converge_on_residual_parameter, c_bool)
    info%converge_on_zero_diff = logical(f_info%converge_on_gradient, c_bool)
    info%fcn_count = f_info%function_evaluation_count
    info%gradient_count = 0
    info%iter_count = f_info%iteration_count
    info%jacobian_count = 0

    do i = 1, nparams
        stats(i)%confidence_interval = f_stats(i)%confidence_interval
        stats(i)%probability = f_stats(i)%probability
        stats(i)%standard_error = f_stats(i)%standard_error
        stats(i)%t_statistic = f_stats(i)%t_statistic
    end do
end subroutine

! --------------------
subroutine siso_fit_ode(t, x, dxdt, args)
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: dxdt
    class(*), intent(inout), optional :: args
    class(*), pointer :: ptr
    real(real64) :: frc
    real(real64), allocatable, dimension(:) :: mdl
    select type (args)
    class is (model_information)
        ptr => args%user_info
        mdl = args%model
        frc = args%excitation%interpolate_value(t)
        select type (ptr)
        class is (c_siso_fit_container)
            call ptr%odefcn(size(x), size(mdl), mdl, t, x, frc, dxdt)
        end select
    end select
end subroutine

! --------------------
subroutine siso_constraint_equations(xg, fg, xc, p, fc, args)
    real(real64), intent(in), dimension(:) :: xg
    real(real64), intent(in), dimension(:) :: fg
    real(real64), intent(in), dimension(:) :: xc
    real(real64), intent(in), dimension(:) :: p
    real(real64), intent(out), dimension(:) :: fc
    class(*), intent(inout), optional :: args

    select type (args)
    class is (c_siso_fit_container)
        call args%constraints(size(xg), size(xc), size(p), xg, fg, xc, p, fc)
    end select
end subroutine

! ------------------------------------------------------------------------------
subroutine c_set_lm_solver_options_defaults(x) &
    bind(C, name = "c_set_lm_solver_options_defaults")
    type(c_lm_solver_options), intent(inout) :: x
    type(lm_solver_options) :: opt
    call opt%set_to_default()
    x%damping_decrease_factor = opt%damping_decrease_factor
    x%damping_increase_factor = opt%damping_increase_factor
    x%finite_difference_step_size = opt%finite_difference_step_size
    x%method = opt%method
end subroutine

! ******************************************************************************
! DYNAMICS_QUATERNIONS
! ------------------------------------------------------------------------------
pure subroutine convert_to_c_quaternion(qc, qf)
    type(c_quaternion), intent(out) :: qc
    type(quaternion), intent(in) :: qf
    qc%w = qf%w
    qc%x = qf%x
    qc%y = qf%y
    qc%z = qf%z
end subroutine

! ------------------------------------------------------------------------------
pure subroutine convert_from_c_quaternion(qf, qc)
    type(quaternion), intent(out) :: qf
    type(c_quaternion), intent(in) :: qc
    qf%w = qc%w
    qf%x = qc%x
    qf%y = qc%y
    qf%z = qc%z
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_from_array(x, q) &
    bind(C, name = "c_quaternion_from_array")
    real(c_double), intent(in) :: x(4)
    type(c_quaternion), intent(out) :: q
    q%w = x(1)
    q%x = x(2)
    q%y = x(3)
    q%z = x(4)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_from_matrix(x, ldx, q) &
    bind(C, name = "c_quaternion_from_matrix")
    integer(c_int), intent(in), value :: ldx
    real(c_double), intent(in) :: x(ldx,3)
    type(c_quaternion), intent(out) :: q

    q = quaternion(x(1:3,1:3))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_from_angle_axis(angle, axis, q) &
    bind(C, name = "c_quaternion_from_angle_axis")
    real(c_double), intent(in), value :: angle
    real(c_double), intent(in) :: axis(3)
    type(c_quaternion), intent(out) :: q

    q = quaternion(angle, axis)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_normalize(q) bind(C, name = "c_quaternion_normalize")
    type(c_quaternion), intent(inout) :: q
    type(quaternion) :: qf
    qf = q
    call qf%normalize()
    q = qf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_add(x, y, q) bind(C, name = "c_quaternion_add")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf + yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_subtract(x, y, q) bind(C, name = "c_quaternion_subtract")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf - yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_multiply(x, y, q) bind(C, name = "c_quaternion_multiply")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf * yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_divide(x, y, q) bind(C, name = "c_quaternion_divide")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf / yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_scale(x, y, q) bind(C, name = "c_quaternion_scale")
    real(c_double), intent(in), value :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: yf
    yf = y
    q = x * yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_conjugate(q, qc) bind(C, name = "c_quaternion_conjugate")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: qc
    type(quaternion) :: qf
    qf = q
    qc = conjg(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_rotate(q, r, rp) bind(C, name = "c_quaternion_rotate")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(in) :: r(3)
    real(c_double), intent(out) :: rp(3)
    type(quaternion) :: qf
    qf = q
    rp = aimag(qf * r * conjg(qf))
end subroutine

! ------------------------------------------------------------------------------
function c_quaternion_abs(q) result(rst) bind(C, name = "c_quaternion_abs")
    type(c_quaternion), intent(in) :: q
    real(c_double) :: rst
    type(quaternion) :: qf
    qf = q
    rst = abs(qf)
end function

! ------------------------------------------------------------------------------
subroutine c_quaternion_inverse(q, qinv) bind(C, name = "c_quaternion_inverse")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: qinv
    type(quaternion) :: qf
    qf = q
    qinv = inverse(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_to_matrix(q, r, ldr) bind(C, name = "c_quaternion_to_matrix") 
    type(c_quaternion), intent(in) :: q
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,3)
    type(quaternion) :: qf
    qf = q
    r(1:3,1:3) = qf%to_matrix()
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_to_angle_axis(q, angle, axis) &
    bind(C, name = "c_quaternion_to_angle_axis")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(out) :: angle
    real(c_double), intent(out) :: axis(3)
    type(quaternion) :: qf
    qf = q
    call qf%to_angle_axis(angle, axis)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_exp(q, rst) bind(C, name = "c_quaternion_exp")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: rst
    type(quaternion) :: qf
    qf = q
    rst = exp(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_log(q, rst)  bind(C, name = "c_quaternion_log")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: rst
    type(quaternion) :: qf
    qf = q
    rst = log(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_pow(q, exponent, rst) bind(C, name = "c_quaternion_pow")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(in), value :: exponent
    type(c_quaternion), intent(out) :: rst
    type(quaternion) :: qf
    qf = q
    rst = qf**exponent
end subroutine

! ------------------------------------------------------------------------------
function c_quaternion_dot_product(x, y) result(rst) &
    bind(C, name = "c_quaternion_dot_product")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    real(c_double) :: rst
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    rst = dot_product(xf, yf)
end function

! ------------------------------------------------------------------------------
subroutine c_quaternion_to_roll_pitch_yaw(q, roll, pitch, yaw) &
    bind(C, name = "c_quaternion_to_roll_pitch_yaw")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(out) :: roll
    real(c_double), intent(out) :: pitch
    real(c_double), intent(out) :: yaw
    type(quaternion) :: qf
    qf = q
    call qf%to_roll_pitch_yaw(roll, pitch, yaw)
end subroutine

! ******************************************************************************
! DYNAMICS_GEOMETRY.F90
! ------------------------------------------------------------------------------
pure subroutine convert_to_c_line(lc, lf)
    type(c_line), intent(out) :: lc
    type(line), intent(in) :: lf
    lc%r0 = lf%r0
    lc%v = lf%v
end subroutine

! ------------------------------------------------------------------------------
pure subroutine convert_from_c_line(lf, lc)
    type(line), intent(out) :: lf
    type(c_line), intent(in) :: lc
    lf%r0 = lc%r0
    lf%v = lc%v
end subroutine

! ------------------------------------------------------------------------------
pure subroutine convert_to_c_plane(pc, pf)
    type(c_plane), intent(out) :: pc
    type(plane), intent(in) :: pf
    pc%a = pf%a
    pc%b = pf%b
    pc%c = pf%c
    pc%d = pf%d
end subroutine

! ------------------------------------------------------------------------------
pure subroutine convert_from_c_plane(pf, pc)
    type(plane), intent(out) :: pf
    type(c_plane), intent(in) :: pc
    pf%a = pc%a
    pf%b = pc%b
    pf%c = pc%c
    pf%d = pc%d
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_normal(pln, nrm) bind(C, name = "c_plane_normal")
    type(c_plane), intent(in) :: pln
    real(c_double), intent(out) :: nrm(3)
    type(plane) :: p
    p = pln
    nrm = plane_normal(p)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_from_3_points(pt1, pt2, pt3, pln) &
    bind(C, name = "c_plane_from_3_points")
    real(c_double), intent(in) :: pt1(3)
    real(c_double), intent(in) :: pt2(3)
    real(c_double), intent(in) :: pt3(3)
    type(c_plane), intent(out) :: pln
    pln = plane(pt1, pt2, pt3)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_from_point_and_normal(pt, nrm, pln) &
    bind(C, name = "c_plane_from_point_and_normal")
    real(c_double), intent(in) :: pt(3)
    real(c_double), intent(in) :: nrm(3)
    type(c_plane), intent(out) :: pln
    pln = plane(pt, nrm)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_from_points(n, pts, ldp, pln) &
    bind(C, name = "c_plane_from_points")
    integer(c_int), intent(in), value :: n, ldp
    real(c_double), intent(in) :: pts(ldp,3)
    type(c_plane), intent(out) :: pln
    pln = plane(pts(1:n,:))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_flip_plane_normal(pln) bind(C, name = "c_flip_plane_normal")
    type(c_plane), intent(inout) :: pln
    type(plane) :: p
    p = pln
    call p%flip_normal()
    pln = p
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_2_points(pt1, pt2, ln) bind(C, name = "c_line_from_2_points")
    real(c_double), intent(in) :: pt1(3)
    real(c_double), intent(in) :: pt2(3)
    type(c_line), intent(out) :: ln
    ln = line(pt1, pt2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_2_planes(p1, p2, ln) bind(C, name = "c_line_from_2_planes")
    type(c_plane), intent(in) :: p1
    type(c_plane), intent(in) :: p2
    type(c_line), intent(out) :: ln
    type(plane) :: pln1, pln2
    pln1 = p1
    pln2 = p2
    ln = line(pln1, pln2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_points(n, pts, ldp, ln) bind(C, name = "c_line_from_points")
    integer(c_int), intent(in), value :: n, ldp
    real(c_double), intent(in) :: pts(ldp, 3)
    type(c_line), intent(out) :: ln
    ln = line(pts(1:n,:))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_evaluate_line_position(ln, t, x) bind(C, name = "c_evaluate_line_position")
    type(c_line), intent(in) :: ln
    real(c_double), intent(in), value :: t
    real(c_double), intent(out) :: x(3)
    type(line) :: l
    l = ln
    x = l%evaluate(t)
end subroutine

! ------------------------------------------------------------------------------
function c_is_parallel_vectors(n, x, y, tol) result(rst) &
    bind(C, name = "c_is_parallel_vectors")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    rst = logical(is_parallel(x, y, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_parallel_lines(x, y, tol) result(rst) &
    bind(C, name = "c_is_parallel_lines")
    type(c_line), intent(in) :: x
    type(c_line), intent(in) :: y
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(line) :: xf, yf
    xf = x
    yf = y
    rst = logical(is_parallel(xf, yf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_parallel_planes(x, y, tol) result(rst) &
    bind(C, name = "c_is_parallel_planes")
    type(c_plane), intent(in) :: x
    type(c_plane), intent(in) :: y
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(plane) :: xf, yf
    xf = x
    yf = y
    rst = logical(is_parallel(xf, yf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_point_on_plane(pt, pln, tol) result(rst) &
    bind(C, name = "c_is_point_on_plane")
    real(c_double), intent(in) :: pt(3)
    type(c_plane), intent(in) :: pln
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(plane) :: pf
    pf = pln
    rst = logical(is_point_on_plane(pt, pf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_point_on_line(pt, ln, tol) result(rst) &
    bind(C, name = "c_is_point_on_line")
    real(c_double), intent(in) :: pt(3)
    type(c_line), intent(in) :: ln
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(line) :: lf
    lf = ln
    rst = logical(is_point_on_line(pt, lf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_nearest_point_on_line(pt, ln) result(rst) &
    bind(C, name = "c_nearest_point_on_line")
    real(c_double), intent(in) :: pt(3)
    type(c_line), intent(in) :: ln
    real(c_double) :: rst
    type(line) :: lf
    lf = ln
    rst = nearest_point_on_line(pt, lf)
end function

! ------------------------------------------------------------------------------
function c_point_to_line_distance(pt, ln) result(rst) &
    bind(C, name = "c_point_to_line_distance")
    real(c_double), intent(in) :: pt(3)
    type(c_line), intent(in) :: ln
    real(c_double) :: rst
    type(line) :: lf
    lf = ln
    rst = point_to_line_distance(pt, lf)
end function

! ------------------------------------------------------------------------------
function c_point_to_plane_distance(pt, pln) result(rst) &
    bind(C, name = "c_point_to_plane_distance")
    real(c_double), intent(in) :: pt(3)
    type(c_plane), intent(in) :: pln
    real(c_double) :: rst
    type(plane) :: pf
    pf = pln
    rst = point_to_plane_distance(pt, pf)
end function

! ------------------------------------------------------------------------------
subroutine c_vector_plane_projection(x, pln, px) &
    bind(C, name = "c_vector_plane_projection")
    real(c_double), intent(in) :: x(3)
    type(c_plane), intent(in) :: pln
    real(c_double), intent(out) :: px(3)
    type(plane) :: pf
    pf = pln
    px = vector_plane_projection(x, pf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_point_plane_projection(pt, pln, ppt) &
    bind(C, name = "c_point_plane_projection")
    real(c_double), intent(in) :: pt(3)
    type(c_plane), intent(in) :: pln
    real(c_double), intent(out) :: ppt(3)
    type(plane) :: pf
    pf = pln
    ppt = point_plane_projection(pt, pf)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine convert_to_c_plucker_line(pc, pf)
    type(c_plucker_line), intent(out) :: pc
    type(plucker_line), intent(in) :: pf
    pc%u = pf%v(1:3)
    pc%m = pf%v(4:6)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine convert_from_c_plucker_line(pf, pc)
    type(plucker_line), intent(out) :: pf
    type(c_plucker_line), intent(in) :: pc
    pf%v(1:3) = pc%u
    pf%v(4:6) = pc%m
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_2_points(pt1, pt2, ln) &
    bind(C, name = "c_plucker_line_from_2_points")
    real(c_double), intent(in) :: pt1(3)
    real(c_double), intent(in) :: pt2(3)
    type(c_plucker_line), intent(out) :: ln
    ln = plucker_line(pt1, pt2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_line(src, ln) &
    bind(C, name = "c_plucker_line_from_line")
    type(c_line), intent(in) :: src
    type(c_plucker_line), intent(out) :: ln
    type(line) :: fsrc
    fsrc = src
    ln = plucker_line(fsrc)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_2_planes(p1, p2, ln) &
    bind(C, name = "c_plucker_line_from_2_planes")
    type(c_plane), intent(in) :: p1
    type(c_plane), intent(in) :: p2
    type(c_plucker_line), intent(out) :: ln
    type(plane) :: pln1, pln2
    pln1 = p1
    pln2 = p2
    ln = plucker_line(pln1, pln2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_array(x, ln) &
    bind(C, name = "c_plucker_line_from_array")
    real(c_double), intent(in) :: x(6)
    type(c_plucker_line), intent(out) :: ln
    ln = plucker_line(x, nrm = .true.)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_mtx_mult(n, x, ldx, ln, y) &
    bind(C, name = "c_plucker_line_mtx_mult")
    integer(c_int), intent(in), value :: n, ldx
    real(c_double), intent(in) :: x(ldx,6)
    type(c_plucker_line), intent(in) :: ln
    real(c_double), intent(out) :: y(n)
    type(plucker_line) :: l
    l = ln
    y = matmul(x(1:n,:), l)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_to_array(ln, x) &
    bind(C, name = "c_plucker_line_to_array")
    type(c_plucker_line), intent(in) :: ln
    real(c_double), intent(out) :: x(6)
    x(1:3) = ln%u
    x(4:6) = ln%m
end subroutine

! ******************************************************************************
! DYNAMICS_MAPS.F90
! ------------------------------------------------------------------------------
subroutine c_poincare_map(n, x, y, z, pln, side, nbuffer, xbuff, ybuff, zbuff, &
    nactual) bind(C, name = "c_poincare_map")
    use dynamics_maps
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(in) :: z(n)
    type(c_plane), intent(in) :: pln
    integer(c_int), intent(in), value :: side
    integer(c_int), intent(in), value :: nbuffer
    real(c_double), intent(out) :: xbuff(nbuffer)
    real(c_double), intent(out) :: ybuff(nbuffer)
    real(c_double), intent(out) :: zbuff(nbuffer)
    integer(c_int), intent(out) :: nactual

    ! Local Variables
    integer(int32) :: i
    type(plane) :: p
    real(real64), allocatable, dimension(:,:) :: rst

    ! Process
    p = pln
    rst = poincare_map(x, y, z, p, side)
    nactual = min(nbuffer, size(rst, 1))
    do concurrent (i = 1:nactual)
        xbuff(i) = rst(i,1)
        ybuff(i) = rst(i,2)
        zbuff(i) = rst(i,3)
    end do
end subroutine

! ------------------------------------------------------------------------------
end module