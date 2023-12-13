program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use example_models
    implicit none

    ! Parameters
    real(real64), parameter :: fs = 2.56d2
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    type(forced_ode_container) :: mdl
    real(real64), allocatable :: mag(:), phase(:), amag(:), &
        aphase(:)
    complex(real64), allocatable :: tf(:), s(:)
    type(frf) :: rsp

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Define the model
    mdl%fcn => example_2nd_order
    mdl%forcing_function => example_2nd_order_forcing

    ! Compute the FRF's
    rsp = frequency_response(mdl, 5.0d0, [0.0d0, 0.0d0])

    ! Extract the magnitude and phase
    mag = 2.0d1 * log10(abs(rsp%responses(:,1) / rsp%responses(1,1)))
    phase = (1.8d2 / pi) * atan2(aimag(rsp%responses(:,1)), &
        real(rsp%responses(:,1)))

    ! Compute the analytical solution
    s = j * (2.0d0 * pi * rsp%frequency)
    tf = 1.0d0 / (s**2 + 2.0d0 * z * wn * s + wn**2)
    amag = 2.0d1 * log10(abs(tf / tf(1)))
    aphase = (1.8d2 / pi) * atan2(aimag(tf), real(tf))

    ! Plot the results
    call plt%initialize(2, 1)
    call plt1%initialize()
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    lgnd => plt1%get_legend()
    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("|k X / F| [dB]")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 1.0d2)
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(rsp%frequency, mag)
    call pd1%set_line_width(2.0)
    call pd1%set_name("Numerical")
    call plt1%push(pd1)

    call pd2%define_data(rsp%frequency, amag)
    call pd2%set_line_width(2.0)
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_frequency(20)
    call pd2%set_line_style(LINE_DOTTED)
    call pd2%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd2%set_name("Analytical")
    call plt1%push(pd2)

    call plt2%initialize()
    xAxis => plt2%get_x_axis()
    yAxis => plt2%get_y_axis()
    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("{/Symbol f} [rad]")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 1.0d2)

    call pd1%define_data(rsp%frequency, phase)
    call plt2%push(pd1)

    call pd2%define_data(rsp%frequency, aphase)
    call plt2%push(pd2)

    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end program