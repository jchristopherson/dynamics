module dynamics_c_system_id
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
    type(kennedy_carpenter_4), target :: kc4
    type(kennedy_carpenter_5), target :: kc5
    type(tsitouras_54), target :: t54
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
        if (nweights /= nw) error stop DYN_INVALID_INPUT_ERROR
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
    case (DYN_KENNEDY_CARPENTER_4)
        integrator_obj => kc4
    case (DYN_KENNEDY_CARPENTER_5)
        integrator_obj => kc5
    case (DYN_TSITOURAS_5)
        integrator_obj => t54
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

end module
