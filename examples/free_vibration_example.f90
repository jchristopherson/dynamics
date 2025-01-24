program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    implicit none

    ! Model Parameters
    real(real64), parameter :: fn = 1.0d1
    real(real64), parameter :: zeta = 1.0d-1

    ! Initial Conditions
    real(real64), parameter :: xo = 0.0d0
    real(real64), parameter :: vo = 1.0d0

    ! Additional Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    integer(int32), parameter :: n = 1000
    real(real64), parameter :: dt = 1.0d-3

    ! Variables
    integer(int32) :: i
    real(real64) :: t(n), x(n), wn, wd, A, B, fnx, delta, x1, x2, t1, t2, zx

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Generate the response of an SDOF undamped free response
    wn = 2.0d0 * pi * fn
    wd = wn * sqrt(1.0d0 - zeta**2)
    A = (vo + zeta * wn * xo) / wd
    B = xo
    t = (/ (i * dt, i = 0, n - 1) /)
    x = exp(-zeta * wn * t) * (A * sin(wd * t) + B * cos(wd * t))

    ! Attempt to estimate the damping ratio and damped natural frequency
    call find_free_response_properties(t, x, delta, fnx, x1 = x1, x2 = x2, &
        t1 = t1, t2 = t2)
    zx = damping_from_log_decrement(delta)

    ! Write out comparisons
    print "(AF5.3)", "Actual Damping Ratio: ", zeta
    print "(AF5.3)", "Estimated Damping Ratio: ", zx
    print "(AF6.3A)", "Actual Damped Frequency: ", wd / (2.0d0 * pi), " Hz"
    print "(AF6.3A)", "Estimated Damped Frequency: ", fnx, " Hz"

    ! Plot the response
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")

    call pd%define_data(t, x)
    call pd%set_line_width(2.0)
    call plt%push(pd)
    
    call pd%define_data([t1, t2], [x1, x2])
    call pd%set_draw_line(.false.)
    call pd%set_draw_markers(.true.)
    call pd%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd%set_marker_scaling(1.5)
    call plt%push(pd)

    call plt%draw()
end program