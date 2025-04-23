module dynamic_system
    use iso_fortran_env
    use dynamics
    implicit none

    ! System Model:
    ! x" + 2 * z * wn * x' + wn**2 * x = F(t)

    ! Simulation Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: sample_rate = 2.0d2
    real(real64), parameter :: max_time = 1.0d1
    real(real64), parameter :: input_amplitude = 1.0d6
    real(real64), parameter :: min_freq = 1.0d0
    real(real64), parameter :: max_freq = 1.0d2

    ! Model Parameters
    real(real64), parameter :: fn = 5.0d1
    real(real64), parameter :: wn = 2.0d0 * pi * fn
    real(real64), parameter :: zeta = 5.0d-2

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
        dxdt(2) = f - 2.0d0 * zeta * wn * x(2) - wn**2 * x(1)
    end subroutine
end module

program example
    use iso_fortran_env
    use dynamics
    use diffeq
    use fplot_core
    use dynamic_system
    implicit none

    ! Parameters
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Variables
    integer(int32) :: i, n
    real(real64), allocatable, dimension(:) :: t, x, amp, phase, aamp, aphase
    real(real64), allocatable, dimension(:,:) :: y
    complex(real64), allocatable, dimension(:) :: s, arsp
    type(ode_container) :: mdl
    type(runge_kutta_45) :: solver
    type(frf) :: rsp

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: x1, x2, y1, y2
    class(legend), pointer :: lgnd

    ! Define the input signal
    n = floor(max_time * sample_rate) + 1
    allocate(t(n))
    t = (/ (i / sample_rate, i = 0, n - 1) /)
    x = chirp(t, input_amplitude, max_time, min_freq, max_freq)

    ! Compute the response of the dynamic system to the specified input
    mdl%fcn => equations_of_motion
    call solver%solve(mdl, t, [0.0d0, 0.0d0])
    y = solver%get_solution()

    ! Compute the frequency response function
    rsp = frequency_response(x, y(:,2), sample_rate)

    ! Convert from complex form into amplitude & phase
    amp = 2.0d1 * log10(abs(rsp%responses(:,1) / rsp%responses(1,1)))
    phase = 1.8d2 * atan2(aimag(rsp%responses(:,1)), real(rsp%responses(:,1))) / pi

    ! Compute the analytical response for comparison
    s = 2.0d0 * pi * rsp%frequency * j
    arsp = input_amplitude / (s**2 + 2.0d0 * zeta * wn * s + wn**2)
    aamp = 2.0d1 * log10(abs(arsp / arsp(1)))
    aphase = 1.8d2 * atan2(aimag(arsp), real(arsp)) / pi

! ------------------------------------------------------------------------------
! PLOTTING
! ------------------------------------------------------------------------------
    ! Initialize the plot objects
    call plt%initialize(2, 1)
    call plt1%initialize()
    call plt2%initialize()
    x1 => plt1%get_x_axis()
    y1 => plt1%get_y_axis()
    x2 => plt2%get_x_axis()
    y2 => plt2%get_y_axis()
    lgnd => plt1%get_legend()

    ! Define labels
    call x1%set_title("f [Hz]")
    call y1%set_title("Y / X [dB]")
    call x2%set_title("f [Hz]")
    call y2%set_title("{/Symbol f} [deg]")

    ! Set up the legend
    call lgnd%set_is_visible(.true.)

    ! Plot the frequency response function
    call pd1%define_data(rsp%frequency, amp)
    call pd1%set_line_width(2.0)
    call pd1%set_name("SISO")
    call plt1%push(pd1)

    call pd1%define_data(rsp%frequency, phase)
    call plt2%push(pd1)

    ! Plot the analytical response
    call pd2%define_data(rsp%frequency, aamp)
    call pd2%set_line_style(LINE_DASHED)
    call pd2%set_line_width(4.0)
    call pd2%set_name("Analytical")
    call plt1%push(pd2)

    call pd2%define_data(rsp%frequency, aphase)
    call plt2%push(pd2)

    ! Draw the plot
    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end program