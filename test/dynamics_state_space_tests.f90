module dynamics_state_space_tests
    use iso_fortran_env
    use fortran_test_helper
    use diffeq, only : runge_kutta_45, ode_container
    use dynamics
    implicit none

    real(real64), private, parameter :: m = 6.0d-4
    real(real64), private, parameter :: k = 2.0d6
    real(real64), private, parameter :: zeta = 1.0d-1
    real(real64), private, parameter :: b = 2.0d0 * zeta * sqrt(k * m)
    real(real64), private, parameter :: U = 1.0d4
    real(real64), private, parameter :: w = 5.0d2

contains
! ------------------------------------------------------------------------------
subroutine ode_fcn(t, x, dxdt, args)
    real(real64), intent(in) :: t, x(:)
    real(real64), intent(out) :: dxdt(:)
    class(*), intent(inout), optional :: args

    ! m x" + b x' + k x = U sin(w * t)
    dxdt(1) = x(2)
    dxdt(2) = (1.0d0 / m) * (U * sin(w * t) - b * x(2) - k * x(1))
end subroutine

subroutine ode_excite_fcn(t, f, args)
    real(real64), intent(in) :: t
    real(real64), intent(out) :: f(:)
    class(*), intent(inout), optional :: args
    f = U * sin(w * t)
end subroutine

function test_lti_solve() result(rst)
    logical :: rst

    ! Parameters
    integer(int32), parameter :: ntime = 1000
    real(real64), parameter :: tmax = 1.0d-1
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dt, t(ntime), ic(2)
    real(real64), allocatable, dimension(:,:) :: sol, ltiSol
    procedure(ss_excitation), pointer :: excite
    type(ode_container) :: odeMdl
    type(state_space) :: mdl
    type(runge_kutta_45) :: integrator

    ! Initialization
    rst = .true.
    excite => ode_excite_fcn
    odeMdl%fcn => ode_fcn

    ! Build the state space model
    allocate(mdl%A(2, 2), mdl%B(2, 1), mdl%C(1, 2), mdl%D(1, 1), source = 0.0d0)
    mdl%A(2,1) = -k / m
    mdl%A(2,2) = -b / m
    mdl%A(1,2) = 1.0d0
    mdl%B(2,1) = 1.0d0 / m
    mdl%C(1,1) = 1.0d0

    ! Build a time vector and an initial conditions vector
    dt = tmax / (ntime - 1.0d0)
    t = (/ (i * dt, i = 0, ntime - 1) /)
    ic = 0.0d0

    ! Solve using LTI
    ltiSol = lti_solve(mdl, excite, t, ic)

    ! Solve using a traditional ODE solver
    call integrator%solve(odeMdl, t, ic)
    sol = integrator%get_solution()

    ! Test
    if (.not.assert(ltiSol, sol(:,1:2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_lti_solve -1"
    end if
end function

! ------------------------------------------------------------------------------
end module