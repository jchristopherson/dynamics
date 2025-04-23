module duffing_ode_container
    use iso_fortran_env
    use dynamics
    implicit none

    ! Duffing Model Parameters
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 4.0d-2
    real(real64), parameter :: delta = 1.0d-1
    real(real64), parameter :: gamma = 1.0d0

contains
    pure subroutine duffing_ode(freq, x, y, dydx, args)
        real(real64), intent(in) :: freq
            ! The excitation frequency
        real(real64), intent(in) :: x
            ! The independent variable.
        real(real64), intent(in), dimension(:) :: y
            ! An array of the N dependent variables.
        real(real64), intent(out), dimension(:) :: dydx
            ! An output array of length N where the derivatives are written.
        class(*), intent(inout), optional :: args
            ! An optional object for input/output of additional information.

        ! Variables
        real(real64) :: f

        ! Compute the harmonic forcing function
        f = gamma * sin(freq * x)

        ! Compute the derivatives
        dydx(1) = y(2)
        dydx(2) = f - delta * y(2) - alpha * y(1) - beta * y(1)**3
    end subroutine

    ! Analytical Solution - Not necessary for the analysis, but used for 
    ! the purpose of comparison.
    pure elemental function leg1(z) result(rst)
        use ieee_arithmetic
        real(real64), intent(in) :: z
        real(real64) :: rst

        ! Local Variables
        real(real64) :: arg, s, nan

        ! Process
        nan = ieee_value(nan, ieee_quiet_nan)
        arg = 4.0d0 * gamma**2 - 3.0d0 * beta * delta**2 * &
            z**4 + (delta**2 - 4.0d0 * alpha) * delta**2 * z**2
        if (arg < 0.0d0) then
            s = -1.0d0
        else
            s = (2.0d0 * sqrt(arg) + z * (3.0d0 * beta * z**2 - &
                2.0d0 * delta**2 + 4.0d0 * alpha)) / (4.0d0 * z)
        end if
        if (s < 0.0d0) then
            rst = nan
        else
            rst = sqrt(s)
        end if
    end function

    pure elemental function leg2(z) result(rst)
        use ieee_arithmetic
        real(real64), intent(in) :: z
        real(real64) :: rst

        ! Local Variables
        real(real64) :: arg, s, nan

        ! Process
        nan = ieee_value(nan, ieee_quiet_nan)
        arg = 4.0d0 * gamma**2 - 3.0d0 * beta * delta**2 * &
            z**4 + (delta**2 - 4.0d0 * alpha) * delta**2 * z**2
        if (arg < 0.0d0) then
            s = -1.0d0
        else
            s = (-2.0d0 * sqrt(arg) + z * (3.0d0 * beta * z**2 - &
                2.0d0 * delta**2 + 4.0d0 * alpha)) / (4.0d0 * z)
        end if
        if (s < 0.0d0) then
            rst = nan
        else
            rst = sqrt(s)
        end if
    end function
end module


program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use duffing_ode_container
    implicit none

    ! Parameters
    real(real64), parameter :: f1 = 0.5d0
    real(real64), parameter :: f2 = 2.0d0
    integer(int32), parameter :: nfreq = 100
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: zmax = 1.0d1
    real(real64), parameter :: zmin = 0.5d0
    
    ! Local Variables
    procedure(harmonic_ode), pointer :: fcn
    integer(int32) :: i
    real(real64) :: dz, phase1(nfreq), phase2(nfreq), z(npts), w1(npts), w2(npts)
    type(frf) :: solup, soldown

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd1, pd2, pd3, pd4
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Point to the ODE routin
    fcn => duffing_ode

    ! Perform the ascending sweep
    solup = frequency_sweep(fcn, nfreq, f1, f2, [0.0d0, 0.0d0])
    phase1 = atan2(aimag(solup%responses(:,1)), real(solup%responses(:,1)))

    ! Perform the descending sweep
    soldown = frequency_sweep(fcn, nfreq, f2, f1, [0.0d0, 0.0d0])
    phase2 = atan2(aimag(soldown%responses(:,1)), real(soldown%responses(:,1)))

    ! Compute the analytical solution
    dz = (zmax - zmin) / (npts - 1.0d0)
    z = (/ (dz * i + zmin, i = 0, npts - 1) /)
    w1 = leg1(z)
    w2 = leg2(z)

    ! Plot the results
    call plt%initialize(2, 1)
    call plt1%initialize()
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    lgnd => plt1%get_legend()
    call xAxis%set_title("{/Symbol w}")
    call yAxis%set_title("X({/Symbol w})")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_horizontal_position(LEGEND_LEFT)
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 2.0d0)

    call pd1%define_data(solup%frequency, abs(solup%responses(:,1)))
    call pd1%set_name("Ascending")
    call pd1%set_draw_markers(.true.)
    call pd1%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd1%set_line_style(LINE_DOTTED)
    call plt1%push(pd1)

    call pd2%define_data(soldown%frequency, abs(soldown%responses(:,1)))
    call pd2%set_name("Descending")
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_style(MARKER_EMPTY_SQUARE)
    call pd2%set_line_style(LINE_DOTTED)
    call plt1%push(pd2)

    call pd3%define_data(w1, z)
    call pd3%set_name("Analytical")
    call pd3%set_line_color(CLR_BLACK)
    call pd3%set_line_width(2.0)
    call plt1%push(pd3)

    call pd4%define_data(w2, z)
    call pd4%set_line_color(CLR_BLACK)
    call pd4%set_line_width(2.0)
    call plt1%push(pd4)

    call plt2%initialize()
    xAxis => plt2%get_x_axis()
    yAxis => plt2%get_y_axis()
    call xAxis%set_title("{/Symbol w}")
    call yAxis%set_title("{/Symbol f}({/Symbol w})")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 2.0d0)

    call pd1%define_data(solup%frequency, phase1)
    call plt2%push(pd1)

    call pd2%define_data(soldown%frequency, phase2)
    call plt2%push(pd2)

    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end program