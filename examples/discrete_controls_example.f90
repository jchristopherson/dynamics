! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
program discrete_controls_example
    use iso_fortran_env, only : real64
    use dynamics, only : state_space, discrete_state_space, &
        discrete_transfer_function, c2d, d2c, DISCRETIZE_ZOH, DISCRETIZE_TUSTIN
    use fplot_core
    implicit none

    ! A unit force acts on a damped mass, and displacement is measured.
    real(real64), parameter :: mass = 1.0d0
    real(real64), parameter :: stiffness = 16.0d0
    real(real64), parameter :: damping = 1.6d0
    real(real64), parameter :: ts = 0.1d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: wn = sqrt(stiffness / mass)
    real(real64), parameter :: zeta = damping / (2.0d0 * mass * wn)
    real(real64), parameter :: wd = wn * sqrt(1.0d0 - zeta**2)
    integer, parameter :: ns = 81
    integer, parameter :: nf = 400
    integer, parameter :: nc = 361

    integer :: i
    real(real64) :: omega
    real(real64) :: theta
    real(real64), dimension(2) :: x_zoh
    real(real64), dimension(2) :: x_tustin
    real(real64), dimension(ns) :: time
    real(real64), dimension(ns) :: y_zoh
    real(real64), dimension(ns) :: y_tustin
    real(real64), dimension(ns) :: y_tf
    real(real64), dimension(ns) :: y_exact
    real(real64), dimension(nf) :: freq
    real(real64), dimension(nf) :: mag_zoh
    real(real64), dimension(nf) :: mag_tustin
    real(real64), dimension(nf) :: mag_cont
    real(real64), dimension(nf) :: phase_zoh
    real(real64), dimension(nf) :: phase_tustin
    real(real64), dimension(nf) :: phase_cont
    real(real64), dimension(nc) :: circle_x
    real(real64), dimension(nc) :: circle_y
    real(real64), allocatable, dimension(:) :: numerator
    real(real64), allocatable, dimension(:) :: denominator
    complex(real64) :: h_zoh
    complex(real64) :: h_tustin
    complex(real64) :: h_cont
    complex(real64), allocatable, dimension(:) :: poles
    complex(real64), allocatable, dimension(:) :: tf_poles
    complex(real64), allocatable, dimension(:) :: zeros
    type(state_space) :: continuous
    type(state_space) :: recovered
    type(discrete_state_space) :: zoh
    type(discrete_state_space) :: tustin
    type(discrete_transfer_function) :: tf_zoh
    type(discrete_transfer_function) :: tf_tustin
    type(discrete_transfer_function) :: tf_coefficients
    type(multiplot) :: figure
    type(plot_2d) :: step_plot
    type(plot_2d) :: pole_plot
    type(plot_2d) :: magnitude_plot
    type(plot_2d) :: phase_plot
    type(plot_data_2d) :: samples
    type(plot_data_2d) :: pole_data
    type(plot_data_2d) :: zero_data
    class(legend), pointer :: key

    continuous = state_space(mass, damping, stiffness)
    continuous%C = reshape([1.0d0, 0.0d0], [1, 2])
    zoh = c2d(continuous, ts, DISCRETIZE_ZOH)
    tustin = c2d(continuous, ts, DISCRETIZE_TUSTIN)
    tf_zoh = zoh%to_transfer_function()
    tf_tustin = tustin%to_transfer_function()

    ! Coefficients are ordered as b0, b1, b2 and a0, a1, a2 in z**(-1).
    numerator = tf_zoh%Y%get_all()
    denominator = tf_zoh%X%get_all()
    tf_coefficients = discrete_transfer_function(numerator, denominator, ts)
    poles = zoh%poles()
    tf_poles = tf_coefficients%poles()
    zeros = tf_coefficients%zeros()

    print "(A, F6.3, A)", "Sample time: ", ts, " s"
    print "(A, *(ES15.6, 1X))", "Numerator (z^-1):   ", numerator
    print "(A, *(ES15.6, 1X))", "Denominator (z^-1): ", denominator
    print "(A)", "ZOH state-space poles (real, imaginary):"
    do i = 1, size(poles)
        print "(2ES15.6)", real(poles(i), real64), aimag(poles(i))
    end do
    print "(A)", "Transfer-function poles (real, imaginary):"
    do i = 1, size(tf_poles)
        print "(2ES15.6)", real(tf_poles(i), real64), aimag(tf_poles(i))
    end do
    print "(A)", "Transfer-function zeros (real, imaginary):"
    do i = 1, size(zeros)
        print "(2ES15.6)", real(zeros(i), real64), aimag(zeros(i))
    end do
    recovered = d2c(zoh, DISCRETIZE_ZOH)
    print "(A, ES12.4)", "ZOH round-trip max matrix error: ", &
        matrix_error(continuous, recovered)
    recovered = d2c(tustin, DISCRETIZE_TUSTIN)
    print "(A, ES12.4)", "Tustin round-trip max matrix error: ", &
        matrix_error(continuous, recovered)

    ! Both models start from zero state; u(k) = 1 for k >= 0.
    x_zoh = 0.0d0
    x_tustin = 0.0d0
    y_tf = 0.0d0
    do i = 1, ns
        time(i) = real(i - 1, real64) * ts
        y_zoh(i) = dot_product(zoh%C(1,:), x_zoh) + zoh%D(1,1)
        y_tustin(i) = dot_product(tustin%C(1,:), x_tustin) + tustin%D(1,1)
        x_zoh = matmul(zoh%A, x_zoh) + zoh%B(:,1)
        x_tustin = matmul(tustin%A, x_tustin) + tustin%B(:,1)

        ! a0*y(k) = b0*u(k) + b1*u(k-1) + b2*u(k-2)
        !             - a1*y(k-1) - a2*y(k-2); past inputs are zero.
        y_tf(i) = numerator(1)
        if (i > 1) y_tf(i) = y_tf(i) + numerator(2) - denominator(2) * y_tf(i-1)
        if (i > 2) y_tf(i) = y_tf(i) + numerator(3) - denominator(3) * y_tf(i-2)
        y_tf(i) = y_tf(i) / denominator(1)
        y_exact(i) = (1.0d0 - exp(-zeta * wn * time(i)) * &
            (cos(wd * time(i)) + zeta * wn / wd * sin(wd * time(i)))) / stiffness
    end do
    print "(A, ES12.4)", "Max ZOH state-space / difference-equation error: ", &
        maxval(abs(y_zoh - y_tf))
    print "(A, ES12.4)", "Max ZOH / exact sampled step error: ", &
        maxval(abs(y_zoh - y_exact))

    ! Evaluate H(exp(i*omega*Ts)) below Nyquist; the discrete type stores
    ! polynomials but does not provide a frequency-response method.
    do i = 1, nf
        freq(i) = 0.49d0 / ts * real(i - 1, real64) / real(nf - 1, real64)
        omega = 2.0d0 * pi * freq(i)
        h_zoh = frequency_response(tf_coefficients, omega)
        h_tustin = frequency_response(tf_tustin, omega)
        h_cont = 1.0d0 / cmplx(stiffness - mass * omega**2, &
            damping * omega, kind = real64)
        mag_zoh(i) = 20.0d0 * log10(abs(h_zoh))
        mag_tustin(i) = 20.0d0 * log10(abs(h_tustin))
        mag_cont(i) = 20.0d0 * log10(abs(h_cont))
        phase_zoh(i) = atan2(aimag(h_zoh), real(h_zoh, real64)) * 180.0d0 / pi
        phase_tustin(i) = atan2(aimag(h_tustin), real(h_tustin, real64)) * 180.0d0 / pi
        phase_cont(i) = atan2(aimag(h_cont), real(h_cont, real64)) * 180.0d0 / pi
        if (i > 1) then
            if (phase_zoh(i) - phase_zoh(i-1) > 180.0d0) phase_zoh(i) = phase_zoh(i) - 360.0d0
            if (phase_zoh(i) - phase_zoh(i-1) < -180.0d0) phase_zoh(i) = phase_zoh(i) + 360.0d0
        end if
    end do
    do i = 1, nc
        theta = 2.0d0 * pi * real(i - 1, real64) / real(nc - 1, real64)
        circle_x(i) = cos(theta)
        circle_y(i) = sin(theta)
    end do

    call figure%initialize(2, 2, width = 1200, height = 800)
    call step_plot%initialize()
    call pole_plot%initialize()
    call magnitude_plot%initialize()
    call phase_plot%initialize()

    call step_plot%set_title("Unit-step response")
    call step_plot%set_x_axis_title("t [s]")
    call step_plot%set_y_axis_title("Displacement [m]")
    key => step_plot%get_legend()
    call key%set_draw_border(.false.)
    call key%set_is_visible(.true.)
    call step_plot%push(time, y_zoh, name = "ZOH state space", lw = 2.0)
    call step_plot%push(time, y_tustin, name = "Tustin state space", ls = LINE_DASHED)
    call samples%define_data(time, y_tf)
    call samples%set_name("ZOH transfer function")
    call samples%set_draw_line(.false.)
    call samples%set_draw_markers(.true.)
    call samples%set_marker_style(MARKER_EMPTY_CIRCLE)
    call samples%set_marker_scaling(0.6)
    call step_plot%push(samples)
    call figure%set(1, 1, step_plot)

    call pole_plot%set_title("ZOH poles and zeros")
    call pole_plot%set_x_axis_title("Real(z)")
    call pole_plot%set_y_axis_title("Imag(z)")
    call pole_plot%set_square_axes(.true.)
    key => pole_plot%get_legend()
    call key%set_is_visible(.true.)
    call key%set_draw_border(.false.)
    call key%set_horizontal_position(LEGEND_CENTER)
    call key%set_vertical_position(LEGEND_CENTER)
    call pole_plot%push(circle_x, circle_y, name = "Unit circle", ls = LINE_DASHED)
    call pole_data%define_data(real(poles, real64), aimag(poles))
    call pole_data%set_name("Poles")
    call pole_data%set_draw_line(.false.)
    call pole_data%set_draw_markers(.true.)
    call pole_data%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pole_data%set_marker_scaling(1.5)
    call pole_plot%push(pole_data)
    call zero_data%define_data(real(zeros, real64), aimag(zeros))
    call zero_data%set_name("Zeros")
    call zero_data%set_draw_line(.false.)
    call zero_data%set_draw_markers(.true.)
    call zero_data%set_marker_style(MARKER_EMPTY_SQUARE)
    call zero_data%set_marker_scaling(1.5)
    call pole_plot%push(zero_data)
    call figure%set(2, 2, pole_plot)

    call magnitude_plot%set_title("Frequency response magnitude")
    call magnitude_plot%set_x_axis_title("Frequency [Hz]")
    call magnitude_plot%set_y_axis_title("Magnitude [dB]")
    key => magnitude_plot%get_legend()
    call key%set_is_visible(.true.)
    call key%set_draw_border(.false.)
    call magnitude_plot%push(freq, mag_zoh, name = "ZOH")
    call magnitude_plot%push(freq, mag_tustin, name = "Tustin", ls = LINE_DASHED)
    call magnitude_plot%push(freq, mag_cont, name = "Continuous")
    call figure%set(2, 1, magnitude_plot)

    call phase_plot%set_title("Frequency response phase")
    call phase_plot%set_x_axis_title("Frequency [Hz]")
    call phase_plot%set_y_axis_title("Phase [deg]")
    key => phase_plot%get_legend()
    call key%set_is_visible(.true.)
    call key%set_draw_border(.false.)
    call phase_plot%push(freq, phase_zoh, name = "ZOH")
    call phase_plot%push(freq, phase_tustin, name = "Tustin", ls = LINE_DASHED)
    call phase_plot%push(freq, phase_cont, name = "Continuous")
    call figure%set(1, 2, phase_plot)
    call figure%draw()

contains
    function matrix_error(a, b) result(err)
        type(state_space), intent(in) :: a
        type(state_space), intent(in) :: b
        real(real64) :: err

        err = max(maxval(abs(a%A - b%A)), maxval(abs(a%B - b%B)), &
            maxval(abs(a%C - b%C)), maxval(abs(a%D - b%D)))
    end function

    function frequency_response(tf, omega) result(h)
        type(discrete_transfer_function), intent(in) :: tf
        real(real64), intent(in) :: omega
        complex(real64) :: h
        complex(real64) :: delay
        complex(real64) :: num
        complex(real64) :: den
        real(real64), allocatable, dimension(:) :: b
        real(real64), allocatable, dimension(:) :: a
        integer :: j

        b = tf%Y%get_all()
        a = tf%X%get_all()
        delay = exp(cmplx(0.0d0, -omega * tf%sample_time, kind = real64))
        num = (0.0d0, 0.0d0)
        den = (0.0d0, 0.0d0)
        do j = 1, size(b)
            num = num + b(j) * delay**(j - 1)
        end do
        do j = 1, size(a)
            den = den + a(j) * delay**(j - 1)
        end do
        h = num / den
    end function
end program
