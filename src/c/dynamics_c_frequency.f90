module dynamics_c_frequency
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

    if (ldm < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldk < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldms < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldr < nfreq) error stop DYN_INVALID_INPUT_ERROR

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

    if (ldm < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldk < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldms < n) error stop DYN_INVALID_INPUT_ERROR

    
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
    if (ldx < n) error stop DYN_INVALID_INPUT_ERROR
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
    type(kennedy_carpenter_4), target :: kc4
    type(kennedy_carpenter_5), target :: kc5
    type(tsitouras_54), target :: t54
    class(ode_integrator), pointer :: integrator_obj

    type(frf) :: frsp

    if (ldr < nfreq) error stop DYN_INVALID_INPUT_ERROR

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
    case (DYN_KENNEDY_CARPENTER_4)
        integrator_obj => kc4
    case (DYN_KENNEDY_CARPENTER_5)
        integrator_obj => kc5
    case (DYN_TSITOURAS_5)
        integrator_obj => t54
    case default
        integrator_obj => rk45
    end select

    frsp = frequency_sweep(odefcn, freq, iv, solver = integrator_obj, &
        args = arg, ncycles = opts%cycle_count, &
        ntransient = opts%transient_cycles, points = opts%points_per_cycle, &
        inHz = logical(opts%frequency_in_hz))
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
    x%cycle_count = 5
    x%transient_cycles = 30
    x%points_per_cycle = 1000
    x%frequency_in_hz = .false.
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
    procedure(c_window_function), pointer :: cfcn
    integer(c_int) :: m

    
    if (mod(winsize, 2) == 0) then
        m = winsize / 2 + 1
    else
        m = (winsize + 1) / 2
    end if
    if (nf /= m) return

    call c_f_procpointer(winfun, cfcn)
    win%size = winsize
    win%fcn => cfcn

    frsp = frequency_response(x, y, fs, win = win, method = method)
    freq = frsp%frequency
    rsp = frsp%responses(:,1)
end subroutine

end module
