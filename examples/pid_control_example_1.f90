module model
    use iso_fortran_env
    implicit none

    real(real64), public, parameter :: m = 6.0d-4
    real(real64), public, parameter :: k = 2.0d3
    real(real64), public, parameter :: zeta = 3.0d-1
    real(real64), public, parameter :: b = 2.0d0 * zeta * sqrt(k * m)
    real(real64), public, parameter :: X = 0.5d0
    real(real64), public, parameter :: f = 5.0d0
    real(real64), public, parameter :: pi = 2.0d0 * acos(0.0d0)

contains
    subroutine command_signal(t, r, args)
        !! Returns the command signal at time 't'
        real(real64), intent(in) :: t
            !! The current simulation time
        real(real64), intent(out) :: r(:)
            !! The command signal.
        class(*), intent(inout), optional :: args
            !! An optional means of passing in arguments from the outside 
            !! world.

        ! We'll use a square wave as the desired output
        r = sign(X, sin(2.0d0 * pi * f * t))
    end subroutine
end module

program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use model
    implicit none

    ! Control Parameters - Arbitrarily Chosen
    real(real64), parameter :: Kp = 5.0d1
    real(real64), parameter :: Ki = 1.5d1
    real(real64), parameter :: Kd = 2.5d-1
    real(real64), parameter :: tau = 1.0d-2

    ! Solution Parameters
    real(real64), parameter :: tmax = 0.5d0
    real(real64), parameter :: init_position = 0.0d0
    real(real64), parameter :: init_velocity = 0.0d0

    ! Local Variables
    procedure(ss_excitation), pointer :: fcn
    integer(int32) :: i, n
    real(real64) :: ic(4), t(2), ri(1)
    real(real64), allocatable, dimension(:) :: r
    real(real64), allocatable, dimension(:,:) :: sol
    type(state_space) :: plant, mdl

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Define a plant model (simple mass-spring-damper model)
    plant = state_space(m, b, k)
    plant%C = reshape([k * X, 0.0d0], [1, 2])

    ! Create the closed-loop model
    mdl = state_space(Kp, Ki, Kd, tau, plant)

    ! Define the initial conditions and solution time vector
    t = [0.0d0, tmax]
    ic = [init_position, init_velocity, 0.0d0, 0.0d0]

    ! Compute the solution
    fcn => command_signal   ! from the model module
    sol = lti_solve(mdl, fcn, t, ic)

    ! Compute the reference signal at each time point - not needed, but useful
    ! for illustration purposes
    n = size(sol, 1)
    allocate(r(n))
    do i = 1, n
        call fcn(sol(i,1), ri)
        r(i) = ri(1)
    end do

    ! Plot
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    
    call pd%define_data(sol(:,1), sol(:,2))
    call pd%set_name("Output")
    call plt%push(pd)

    call pd%define_data(sol(:,1), r)
    call pd%set_name("Command")
    call plt%push(pd)

    call plt%draw()
end program