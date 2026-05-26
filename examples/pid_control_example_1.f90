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
    real(real64), parameter :: Kp = 8.0d-2
    real(real64), parameter :: Ki = 1.0d3
    real(real64), parameter :: Kd = 5.0d-4
    real(real64), parameter :: tau = 1.0d-3

    ! Solution Parameters
    real(real64), parameter :: tmax = 0.5d0
    real(real64), parameter :: init_position = 0.0d0
    real(real64), parameter :: init_velocity = 0.0d0

    ! Local Variables
    procedure(ss_excitation), pointer :: fcn
    integer(int32) :: i, n
    real(real64) :: ic(4), t(2), ri(1)
    real(real64), allocatable, dimension(:) :: r, freq, omega, mag, phase
    real(real64), allocatable, dimension(:,:) :: sol
    complex(real64), allocatable, dimension(:) :: poles, zeros
    complex(real64), allocatable, dimension(:,:,:) :: Z
    type(state_space) :: plant, mdl

    ! Plot Variables
    type(plot_2d) :: plt, plt1, plt2
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis, x1, y1, x2, y2
    type(multiplot) :: bode

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

    ! Compute the poles and zeros of the closed-loop system
    poles = mdl%poles()
    zeros = mdl%zeros()
    print "(A)", "Poles:"
    do i = 1, size(poles)
        print "(A, A, F0.3, A, F0.3, A)", &
            achar(9), "(", real(poles(i)), ", ", aimag(poles(i)), ")"
    end do
    print "(A)", "Zeros:"
    do i = 1, size(zeros)
        print "(A, A, F0.3, A, F0.3, A)", &
            achar(9), "(", real(zeros(i)), ", ", aimag(zeros(i)), ")"
    end do

    ! Compute the transfer functions
    freq = linspace(1.0d0, 1.0d3, 1000)
    omega = 2.0d0 * pi * freq
    Z = mdl%transfer_function(omega)
    mag = 2.0d1 * log10(abs(Z(1,1,:)))  ! convert to dB
    phase = atan2(aimag(Z(1,1,:)), real(Z(1,1,:))) * 1.8d2 / pi

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
    call pd%set_line_width(2.0)
    call pd%set_name("Output")
    call plt%push(pd)

    call pd%define_data(sol(:,1), r)
    call pd%set_name("Command")
    call plt%push(pd)

    call plt%draw()

    ! Bode Plot
    call bode%initialize(2, 1)
    call plt1%initialize()
    call plt2%initialize()
    x1 => plt1%get_x_axis()
    y1 => plt1%get_y_axis()
    x2 => plt2%get_x_axis()
    y2 => plt2%get_y_axis()
    call x1%set_title("f [Hz]")
    call y1%set_title("Y/U [dB]")
    call x2%set_title("f [Hz]")
    call y2%set_title("{/Symbol f} [deg]")

    call pd%define_data(freq, mag)
    call plt1%push(pd)

    call pd%define_data(freq, phase)
    call plt2%push(pd)

    call bode%set(1, 1, plt1)
    call bode%set(2, 1, plt2)
    call bode%draw()
end program