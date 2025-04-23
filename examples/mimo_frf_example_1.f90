module mdof_dynamic_system
    use iso_fortran_env
    use dynamics
    implicit none

    ! MISO System Model:
    !
    ! m x1" + 2 b x1' - b x2' + 2 k x1 - k x2 = F(t)
    ! m x2" + 2 b x2' - b x1' + 2 k x2 - k x1 = 0

    ! Simulation Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: sample_rate = 6.0d2
    real(real64), parameter :: max_time = 1.0d1
    real(real64), parameter :: input_amplitude = 1.0d9
    real(real64), parameter :: min_freq = 1.0d0
    real(real64), parameter :: max_freq = 0.5d0 * sample_rate

    ! Model Parameters
    real(real64), parameter :: m1 = 4.5d0
    real(real64), parameter :: m2 = 2.0d0
    real(real64), parameter :: k = 2.0d6
    real(real64), parameter :: b = 1.0d2

contains
    pure subroutine equations_of_motion(t, x, dxdt, args)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        class(*), intent(inout), optional :: args

        ! Compute the forcing term
        real(real64) :: f
        f = chirp(t, input_amplitude, max_time, min_freq, max_freq)

        ! Equations
        dxdt(1) = x(2)
        dxdt(2) = (f - b * (2.0d0 * x(2) + x(4)) - k * (2.0d0 * x(1) + x(3))) / m1
        dxdt(3) = x(4)
        dxdt(4) = (-b * (2.0d0 * x(4) + x(2)) - k * (2.0d0 * x(3) + x(1))) / m2
    end subroutine
end module

program example
    use iso_fortran_env
    use dynamics
    use diffeq
    use fplot_core
    use mdof_dynamic_system
    implicit none

    ! Variables
    integer(int32) :: i, n
    real(real64), allocatable, dimension(:) :: t, amp11, amp12, phase11, phase12
    real(real64), allocatable, dimension(:,:) :: x, y, z1
    type(ode_container) :: mdl
    type(runge_kutta_45) :: solver
    type(mimo_frf) :: rsp

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2, plt3, plt4
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: x1, x2, x3, x4, y1, y2, y3, y4
    class(terminal), pointer :: term
    
    ! Define the input signal
    n = floor(max_time * sample_rate) + 1
    allocate(t(n), x(n,1))
    t = (/ (i / sample_rate, i = 0, n - 1) /)
    x(:,1) = chirp(t, input_amplitude, max_time, min_freq, max_freq)

    ! Compute the response of the dynamic system to the specified inputs
    mdl%fcn => equations_of_motion
    call solver%solve(mdl, t, [0.0d0, 0.0d0, 0.0d0, 0.0d0])
    z1 = solver%get_solution()

    ! Compute the frequency response function
    allocate(y(n, 2))
    y(:,1) = z1(:,2)
    y(:,2) = z1(:,4)
    rsp = frequency_response(x, y, sample_rate)

    ! Convert to amplitude
    amp11 = 2.0d1 * log10(abs(rsp%responses(:,1,1) / rsp%responses(1,1,1)))
    amp12 = 2.0d1 * log10(abs(rsp%responses(:,2,1) / rsp%responses(1,2,1)))

    ! Compute the phase terms
    phase11 = 1.8d2 * atan2(aimag(rsp%responses(:,1,1)), &
        real(rsp%responses(:,1,1))) / pi
    phase12 = 1.8d2 * atan2(aimag(rsp%responses(:,2,1)), &
        real(rsp%responses(:,2,1))) / pi

! ------------------------------------------------------------------------------
! PLOTTING
! ------------------------------------------------------------------------------
    ! Initialize the plot objects
    call plt%initialize(2, 2)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()
    x1 => plt1%get_x_axis()
    y1 => plt1%get_y_axis()
    x2 => plt2%get_x_axis()
    y2 => plt2%get_y_axis()
    x3 => plt3%get_x_axis()
    y3 => plt3%get_y_axis()
    x4 => plt4%get_x_axis()
    y4 => plt4%get_y_axis()
    term => plt%get_terminal()

    ! Define the window size
    call term%set_window_width(1000)
    call term%set_window_height(600)

    ! Define Labels
    call x1%set_title("Frequency [Hz]")
    call y1%set_title("Magnitude [dB]")
    call x2%set_title("Frequency [Hz]")
    call y2%set_title("Magnitude [dB]")
    call x3%set_title("Frequency [Hz]")
    call y3%set_title("Phase [deg]")
    call x4%set_title("Frequency [Hz]")
    call y4%set_title("Phase [deg]")

    ! Define Titles
    call plt1%set_title("Input 1, Output 1")
    call plt2%set_title("Input 1, Output 2")

    ! Plot
    call pd%define_data(rsp%frequency, amp11)
    call pd%set_line_width(2.0)
    call plt1%push(pd)
    call plt%set(1, 1, plt1)

    call pd%define_data(rsp%frequency, amp12)
    call plt2%push(pd)
    call plt%set(1, 2, plt2)

    call pd%define_data(rsp%frequency, phase11)
    call plt3%push(pd)
    call plt%set(2, 1, plt3)

    call pd%define_data(rsp%frequency, phase12)
    call plt4%push(pd)
    call plt%set(2, 2, plt4)

    call plt%draw()
end program
