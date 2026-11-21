! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
program example
    ! This example illustrates the helper routines in dynamics_vibrations.  A
    ! free response is used to estimate the damping ratio and damped frequency
    ! from peaks and the logarithmic decrement.  A step response is used to
    ! evaluate the rise time, estimate the settling amplitude, and estimate the
    ! damping ratio from fractional overshoot.  Finally, the Q-factor and
    ! bandwidth are evaluated over a range of damping ratios.
    use iso_fortran_env
    use dynamics
    use fplot_core
    implicit none

    ! Model Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: fn = 2.0d0
    real(real64), parameter :: wn = 2.0d0 * pi * fn

    ! Free Response Parameters
    real(real64), parameter :: zeta_free = 5.0d-2
    real(real64), parameter :: dt_free = 1.0d-3
    integer, parameter :: n_free = 5001

    ! Step Response Parameters
    real(real64), parameter :: zeta_step = 3.0d-1
    real(real64), parameter :: xs = 1.0d0
    real(real64), parameter :: dt_step = 1.0d-2
    integer, parameter :: n_step = 2001

    ! Sweep Parameters
    real(real64), parameter :: zeta_min = 5.0d-2
    real(real64), parameter :: zeta_max = 1.0d0
    integer, parameter :: n_sweep = 100

    ! Local Variables
    integer :: i
    real(real64) :: wd, delta, delta_exact, fnx, x1, x2, t1, t2, zeta_est
    real(real64) :: xf, zeta_os, tr, y_tr, dz
    real(real64) :: q_free, bw_free
    real(real64), allocatable, dimension(:) :: t_free
    real(real64), allocatable, dimension(:) :: x_free
    real(real64), allocatable, dimension(:) :: t_step
    real(real64), allocatable, dimension(:) :: y_step
    real(real64), allocatable, dimension(:) :: zeta_sweep
    real(real64), allocatable, dimension(:) :: q_sweep
    real(real64), allocatable, dimension(:) :: bw_sweep

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: p_free, p_step, p_q, p_bw
    type(plot_data_2d) :: pd_peaks, pd_rise, pd_q_marker

! ******************************************************************************
! Free Response: Peaks, Logarithmic Decrement, and Damped Frequency
    wd = wn * sqrt(1.0d0 - zeta_free**2)
    allocate(t_free(n_free), x_free(n_free))
    t_free = [(real(i - 1, real64) * dt_free, i = 1, n_free)]
    x_free = exp(-zeta_free * wn * t_free) * sin(wd * t_free) / wd

    call find_free_response_properties(t_free, x_free, delta, fnx, &
        x1 = x1, x2 = x2, t1 = t1, t2 = t2)
    zeta_est = damping_from_log_decrement(delta)
    delta_exact = 2.0d0 * pi * zeta_free / sqrt(1.0d0 - zeta_free**2)

    print "(A)", "Free response"
    print "(A, F8.4)", "    Actual logarithmic decrement:     ", delta_exact
    print "(A, F8.4)", "    Estimated logarithmic decrement:  ", delta
    print "(A, F8.4)", "    Actual damping ratio:             ", zeta_free
    print "(A, F8.4)", "    Estimated damping ratio:          ", zeta_est
    print "(A, F8.4, A)", "    Actual damped frequency:          ", &
        wd / (2.0d0 * pi), " Hz"
    print "(A, F8.4, A)", "    Estimated damped frequency:       ", fnx, " Hz"
    print "(A)", ""

! ******************************************************************************
! Step Response: Rise Time, Settling Amplitude, and Fractional Overshoot
    allocate(t_step(n_step), y_step(n_step))
    t_step = [(real(i - 1, real64) * dt_step, i = 1, n_step)]
    y_step = evaluate_step_response(wn, zeta_step, xs, t_step)

    xf = find_settling_amplitude(y_step)
    zeta_os = damping_from_fractional_overshoot(y_step)
    tr = rise_time(wn, zeta_step)
    y_tr = evaluate_step_response(wn, zeta_step, xs, tr)

    print "(A)", "Step response"
    print "(A, F8.4, A)", "    Rise time:                        ", tr, " s"
    print "(A, F8.4)", "    Response at rise time:            ", y_tr
    print "(A, F8.4)", "    Settling amplitude:               ", xf
    print "(A, F8.4)", "    Actual damping ratio:             ", zeta_step
    print "(A, F8.4)", "    Estimated damping ratio (overshoot): ", zeta_os
    print "(A)", ""

! ******************************************************************************
! Q-Factor and Bandwidth Over a Range of Damping Ratios
    dz = (zeta_max - zeta_min) / real(n_sweep - 1, real64)
    allocate(zeta_sweep(n_sweep), q_sweep(n_sweep), bw_sweep(n_sweep))
    zeta_sweep = [(zeta_min + real(i - 1, real64) * dz, i = 1, n_sweep)]
    q_sweep = q_factor(zeta_sweep)
    bw_sweep = estimate_bandwidth(fn, zeta_sweep)

    q_free = q_factor(zeta_free)
    bw_free = estimate_bandwidth(fn, zeta_free)
    print "(A)", "Q-factor and bandwidth at the free-response damping ratio"
    print "(A, F8.4)", "    Q-factor:                         ", q_free
    print "(A, F8.4, A)", "    Bandwidth:                        ", bw_free, " Hz"

! ******************************************************************************
! Plotting
    call plt%initialize(2, 2, width = 1200, height = 800)
    call p_free%initialize()
    call p_step%initialize()
    call p_q%initialize()
    call p_bw%initialize()

    ! Free response with detected peaks
    call p_free%set_title("Free response")
    call p_free%set_x_axis_title("t [s]")
    call p_free%set_y_axis_title("x(t)")
    call p_free%push(t_free, x_free, lw = 2.0)
    call pd_peaks%define_data([t1, t2], [x1, x2])
    call pd_peaks%set_draw_line(.false.)
    call pd_peaks%set_draw_markers(.true.)
    call pd_peaks%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd_peaks%set_marker_scaling(1.5)
    call p_free%push(pd_peaks)
    call plt%set(1, 1, p_free)

    ! Step response with settling amplitude and rise time
    call p_step%set_title("Step response")
    call p_step%set_x_axis_title("t [s]")
    call p_step%set_y_axis_title("x(t) / X_s")
    call p_step%push(t_step, y_step, lw = 2.0)
    call p_step%push([t_step(1), t_step(n_step)], [xf, xf], ls = LINE_DASHED)
    call pd_rise%define_data([tr], [y_tr])
    call pd_rise%set_draw_line(.false.)
    call pd_rise%set_draw_markers(.true.)
    call pd_rise%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd_rise%set_marker_scaling(1.5)
    call p_step%push(pd_rise)
    call plt%set(2, 1, p_step)

    ! Q-factor versus damping ratio
    call p_q%set_title("Q-factor")
    call p_q%set_x_axis_title("{/Symbol z}")
    call p_q%set_y_axis_title("Q")
    call p_q%push(zeta_sweep, q_sweep, lw = 2.0)
    call pd_q_marker%define_data([zeta_free], [q_free])
    call pd_q_marker%set_draw_line(.false.)
    call pd_q_marker%set_draw_markers(.true.)
    call pd_q_marker%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd_q_marker%set_marker_scaling(1.5)
    call p_q%push(pd_q_marker)
    call plt%set(1, 2, p_q)

    ! Bandwidth versus damping ratio
    call p_bw%set_title("Bandwidth (f_n = 2 Hz)")
    call p_bw%set_x_axis_title("{/Symbol z}")
    call p_bw%set_y_axis_title("{/Symbol D}f [Hz]")
    call p_bw%push(zeta_sweep, bw_sweep, lw = 2.0)
    call plt%set(2, 2, p_bw)

    call plt%draw()
end program