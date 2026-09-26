module dynamics_c_controls
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

! ------------------------------------------------------------------------------
subroutine c_evaluate_transfer_function(tf, n, s, z) &
    bind(C, name = "c_evaluate_transfer_function")
    type(c_transfer_function), intent(in) :: tf
    integer(c_int), intent(in), value :: n
    complex(c_double), intent(in) :: s(n)
    complex(c_double), intent(out) :: z(n)

    type(transfer_function) :: t
    t = tf
    z = t%evaluate(s)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_transfer_function_poles(tf, n, p) &
    bind(C, name = "c_transfer_function_poles")
    type(c_transfer_function), intent(in) :: tf
    integer(c_int), intent(in), value :: n
    complex(c_double), intent(out) :: p(n)

    integer(int32) :: np, mnp
    complex(real64), allocatable, dimension(:) :: poles
    type(transfer_function) :: t
    t = tf
    poles = t%poles()
    np = size(poles)
    mnp = min(n, np)
    p(1:mnp) = poles(1:mnp)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_transfer_function_zeros(tf, n, z) &
    bind(C, name = "c_transfer_function_zeros")
    type(c_transfer_function), intent(in) :: tf
    integer(c_int), intent(in), value :: n
    complex(c_double), intent(out) :: z(n)

    integer(int32) :: nz, mnz
    complex(real64), allocatable, dimension(:) :: zeros
    type(transfer_function) :: t
    t = tf
    zeros = t%zeros()
    nz = size(zeros)
    mnz = min(nz, n)
    z(1:mnz) = zeros(1:mnz)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_to_ccf_state_space(tf, ss) bind(C, name = "c_to_ccf_state_space")
    type(c_transfer_function), intent(in) :: tf
    type(c_state_space_model), intent(out) :: ss

    type(transfer_function) :: t
    type(state_space) :: s
    t = tf
    s = t%to_ccf_state_space()
    ss = s
end subroutine

! ------------------------------------------------------------------------------
subroutine c_to_ocf_state_space(tf, ss) bind(C, name = "c_to_ocf_state_space")
    type(c_transfer_function), intent(in) :: tf
    type(c_state_space_model), intent(out) :: ss

    type(transfer_function) :: t
    type(state_space) :: s
    t = tf
    s = t%to_ocf_state_space()
    ss = s
end subroutine

! ------------------------------------------------------------------------------
subroutine c_create_state_space_model(n, n_out, m, ldm, b, ldb, k, ldk, mdl) &
    bind(C, name = "c_create_state_space_model")
    integer(c_int), intent(in), value :: n, n_out, ldm, ldb, ldk
    real(c_double), intent(in) :: m(ldm,n), b(ldb,n), k(ldk,n)
    type(c_state_space_model), intent(out) :: mdl
    type(state_space) :: ss
    if (ldm < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldb < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldk < n) error stop DYN_INVALID_INPUT_ERROR
    ss = state_space(m(1:n,:), b(1:n,:), k(1:n,:), n_out)
    mdl = ss
end subroutine

! ------------------------------------------------------------------------------
subroutine c_create_pid_state_space_model(kp, ki, kd, tau, plant, mdl) &
    bind(C, name = "c_create_pid_state_space_model")
    real(c_double), intent(in), value :: kp, ki, kd, tau
    type(c_state_space_model), intent(in) :: plant
    type(c_state_space_model), intent(out) :: mdl

    type(state_space) :: fplant, fmdl
    fplant = plant
    fmdl = state_space(kp, ki, kd, tau, fplant)
    mdl = fmdl
end subroutine

! ------------------------------------------------------------------------------
subroutine c_transfer_function_multiply(tf1, tf2, tf) &
    bind(C, name = "c_transfer_function_multiply")
    type(c_transfer_function), intent(in) :: tf1, tf2
    type(c_transfer_function), intent(out) :: tf
    type(transfer_function) :: t1, t2, t
    t1 = tf1
    t2 = tf2
    t = t1 * t2
    tf = t
end subroutine

! ------------------------------------------------------------------------------
subroutine c_scale_transfer_function(x, tf1, tf) &
    bind(C, name = "c_scale_transfer_function")
    real(c_double), intent(in), value :: x
    type(c_transfer_function), intent(in) :: tf1
    type(c_transfer_function), intent(out) :: tf
    type(transfer_function) :: t1, t
    t1 = tf1
    t = x * t1
    tf = t
end subroutine

! ------------------------------------------------------------------------------
subroutine c_lti_solve(mdl, u, n, t, ndof, ic, solver, nout, y, ldy) &
    bind(C, name = "c_lti_solve")
    type(c_state_space_model), intent(in) :: mdl
    type(c_funptr), intent(in), value :: u
    integer(c_int), intent(in), value :: n, ndof, ldy, solver, nout
    real(c_double), intent(in) :: t(n)
    real(c_double), intent(in) :: ic(ndof)
    real(c_double), intent(out) :: y(ldy,nout)

    type(c_ss_excitation_container) :: arg
    procedure(c_ss_excitation), pointer :: fptr
    procedure(ss_excitation), pointer :: ptr
    real(real64), allocatable, dimension(:,:) :: sol
    class(ode_integrator), pointer :: integrator
    type(runge_kutta_23), target :: rk23
    type(runge_kutta_45), target :: rk45
    type(runge_kutta_853), target :: rk853
    type(rosenbrock), target :: rbk
    type(adams), target :: adms
    type(bdf), target :: bdiff
    type(kennedy_carpenter_4), target :: kc4
    type(kennedy_carpenter_5), target :: kc5
    type(tsitouras_54), target :: t54
    type(state_space) :: fmdl

    if (n <= 2) error stop DYN_INVALID_INPUT_ERROR
    if (ldy < n) error stop DYN_INVALID_INPUT_ERROR

    call c_f_procpointer(u, fptr)
    arg%fcn => fptr
    ptr => c_lti_solver_routine

    fmdl = mdl

    select case (solver)
    case (DYN_RUNGE_KUTTA_23)
        integrator => rk23
    case (DYN_RUNGE_KUTTA_45)
        integrator => rk45
    case (DYN_RUNGE_KUTTA_853)
        integrator => rk853
    case (DYN_ROSENBROCK)
        integrator => rbk
    case (DYN_BDF)
        integrator => bdiff
    case (DYN_ADAMS)
        integrator => adms
    case (DYN_KENNEDY_CARPENTER_4)
        integrator => kc4
    case (DYN_KENNEDY_CARPENTER_5)
        integrator => kc5
    case (DYN_TSITOURAS_5)
        integrator => t54
    case default
        integrator => rk45
    end select

    sol = lti_solve(fmdl, ptr, t, ic, solver = integrator, args = arg)
    y(1:n,1:nout) = sol(1:n,2:nout+1)
end subroutine

subroutine c_lti_solver_routine(t, u, args)
    real(real64), intent(in) :: t
    real(real64), intent(out), dimension(:) :: u
    class(*), intent(inout), optional :: args

    integer(int32) :: n
    n = size(u)
    select type (args)
    class is (c_ss_excitation_container)
        call args%fcn(n, t, u)
    end select
end subroutine

! ------------------------------------------------------------------------------
subroutine c_state_space_poles(mdl, n, p) bind(C, name = "c_state_space_poles")
    type(c_state_space_model), intent(in) :: mdl
    integer(c_int), intent(in), value :: n
    complex(c_double), intent(out) :: p(n)

    integer(c_int) :: np
    complex(real64), allocatable, dimension(:) :: poles
    type(state_space) :: fmdl

    fmdl = mdl
    poles = fmdl%poles()
    np = min(size(poles), n)
    p(1:np) = poles(1:np)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_state_space_zeros(mdl, n, z, nz) &
    bind(C, name = "c_state_space_zeros")
    type(c_state_space_model), intent(in) :: mdl
    integer(c_int), intent(in), value :: n
    complex(c_double), intent(out) :: z(n)
    integer(c_int), intent(out) :: nz

    complex(real64), allocatable, dimension(:) :: zeros
    type(state_space) :: fmdl

    fmdl = mdl
    zeros = fmdl%zeros()
    nz = min(size(zeros), n)
    z(1:nz) = zeros(1:nz)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_state_space_transfer_function(mdl, nin, nout, n, s, z, ldz) &
    bind(C, name = "c_state_space_transfer_function")
    type(c_state_space_model), intent(in) :: mdl
    integer(c_int), intent(in), value :: nin, nout, n, ldz
    complex(c_double), intent(in) :: s(n)
    complex(c_double), intent(out) :: z(ldz,nout,n)

    type(state_space) :: fmdl

    fmdl = mdl

    if (nin /= size(fmdl%B, 2)) error stop DYN_INVALID_INPUT_ERROR
    if (nout /= size(fmdl%C, 1)) error stop DYN_INVALID_INPUT_ERROR
    if (ldz < nin) error stop DYN_INVALID_INPUT_ERROR

    z(1:nin,:,:) = fmdl%transfer_function(s)
end subroutine

end module
