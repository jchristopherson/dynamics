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
module mdof_dynamic_system
    use iso_fortran_env, only : real64
    implicit none

    ! Each mass is attached to ground and to the other mass by a spring
    ! and damper. The input and output channels are ordered by mass.
    real(real64), parameter :: sample_rate = 2048.0d0
    real(real64), parameter :: m1 = 4.5d0
    real(real64), parameter :: m2 = 2.0d0
    real(real64), parameter :: k = 2.0d6
    real(real64), parameter :: b = 1.0d2
    real(real64), parameter, dimension(2,2) :: mass = &
        reshape([m1, 0.0d0, 0.0d0, m2], [2,2])
    real(real64), parameter, dimension(2,2) :: stiffness = &
        k * reshape([2.0d0, -1.0d0, -1.0d0, 2.0d0], [2,2])
    real(real64), parameter, dimension(2,2) :: damping = (b / k) * stiffness

    type force_history
        real(real64), allocatable, dimension(:,:) :: samples
    end type

contains
    subroutine equations_of_motion(t, x, dxdt, args)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        class(*), intent(inout), optional :: args
        real(real64), dimension(2) :: force
        real(real64) :: sample
        real(real64) :: fraction
        integer :: i

        ! Interpolate the same force samples used by the FRF estimator.
        if (.not.present(args)) error stop "Missing force history."
        select type (args)
        type is (force_history)
            sample = min(max(t * sample_rate, 0.0d0), &
                real(size(args%samples, 1) - 1, real64))
            i = min(int(sample) + 1, size(args%samples, 1) - 1)
            fraction = sample - real(i - 1, real64)
            force = (1.0d0 - fraction) * args%samples(i,:) + &
                fraction * args%samples(i+1,:)
        class default
            error stop "Expected a force history."
        end select

        ! State: [displacement 1, displacement 2, velocity 1, velocity 2].
        dxdt(1:2) = x(3:4)
        dxdt(3:4) = (force - matmul(damping, x(3:4)) - &
            matmul(stiffness, x(1:2))) / [m1, m2]
    end subroutine
end module

program example
    use iso_fortran_env, only : real64
    use dynamics, only : dynamic_stiffness, frequency_response, mimo_frf
    use diffeq, only : ode_container, runge_kutta_45
    use spectrum, only : hann_window
    use fplot_core
    use mdof_dynamic_system
    implicit none

    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: force_amplitude = 1.0d3
    real(real64), parameter :: min_freq = 10.0d0
    real(real64), parameter :: max_freq = 300.0d0
    integer, parameter :: n = 262144
    integer, parameter :: window_size = 4096

    integer :: i
    integer :: input
    integer :: output
    integer :: panel
    integer :: seed_size
    integer, allocatable, dimension(:) :: seed
    real(real64), allocatable, dimension(:) :: t
    real(real64), allocatable, dimension(:,:) :: solution
    real(real64), allocatable, dimension(:) :: frequency
    complex(real64), dimension(2,2) :: impedance
    complex(real64) :: determinant
    complex(real64), allocatable, dimension(:,:,:) :: exact
    complex(real64), allocatable, dimension(:) :: estimate
    complex(real64), allocatable, dimension(:) :: reference
    character(len = 32) :: title
    type(force_history) :: forces
    type(ode_container) :: model
    type(runge_kutta_45) :: solver
    type(hann_window) :: win
    type(mimo_frf) :: rsp
    type(plot_2d), dimension(2,4) :: axes
    type(multiplot) :: plt

    ! Independent broadband forces let the H1 estimate distinguish the
    ! channels. Window averaging reduces finite-record cross-correlation.
    call random_seed(size = seed_size)
    allocate(seed(seed_size), source = 173)
    call random_seed(put = seed)
    allocate(forces%samples(n,2), t(n))
    call random_number(forces%samples)
    forces%samples = force_amplitude * (2.0d0 * forces%samples - 1.0d0)
    do i = 1, n
        t(i) = (i - 1) / sample_rate
    end do

    model%fcn => equations_of_motion
    ! The adaptive solver may take several steps per force sample.
    call solver%set_step_limit(20 * n)
    call solver%solve(model, t, [0.0d0, 0.0d0, 0.0d0, 0.0d0], args = forces)
    solution = solver%get_solution()

    win%size = window_size
    rsp = frequency_response(forces%samples, solution(:,2:3), sample_rate, win = win)

    ! Compare the estimates around the modes with H = (K - w^2 M + i w C)^-1.
    frequency = pack(rsp%frequency, in_band(rsp%frequency))
    allocate(exact(2,2,size(frequency)))
    do i = 1, size(frequency)
        call dynamic_stiffness(2.0d0 * pi * frequency(i), mass, damping, stiffness, impedance)
        determinant = impedance(1,1) * impedance(2,2) - impedance(1,2) * impedance(2,1)
        exact(:,:,i) = reshape([impedance(2,2), -impedance(2,1), &
            -impedance(1,2), impedance(1,1)], [2,2]) / determinant
    end do

    call plt%initialize(2, 4, width = 1400, height = 650)
    do input = 1, 2
        do output = 1, 2
            panel = 2 * (input - 1) + output
            estimate = pack(rsp%responses(output,input,:), in_band(rsp%frequency))
            reference = exact(output,input,:)
            write(title, '("Input ", I1, " -> Output ", I1)') input, output
            print '(A, ": relative RMS error = ", F5.1, "%")', trim(title), &
                100.0d0 * norm2(abs(estimate - reference)) / norm2(abs(reference))
            call add_panel(1, panel, "Magnitude [dB re 1 m/N]", &
                magnitude_db(estimate), magnitude_db(reference), title)
            call add_panel(2, panel, "Phase [deg]", phase_deg(estimate), phase_deg(reference))
        end do
    end do
    call plt%draw()

contains
    elemental logical function in_band(frequency)
        real(real64), intent(in) :: frequency
        in_band = frequency >= min_freq .and. frequency <= max_freq
    end function

    elemental real(real64) function magnitude_db(h)
        complex(real64), intent(in) :: h
        magnitude_db = 20.0d0 * log10(max(abs(h), tiny(1.0d0)))
    end function

    elemental real(real64) function phase_deg(h)
        complex(real64), intent(in) :: h
        phase_deg = atan2(aimag(h), real(h, real64)) * 180.0d0 / pi
    end function

    subroutine add_panel(row, col, label, estimated_values, exact_values, panel_title)
        integer, intent(in) :: row
        integer, intent(in) :: col
        character(len = *), intent(in) :: label
        real(real64), intent(in), dimension(:) :: estimated_values
        real(real64), intent(in), dimension(:) :: exact_values
        character(len = *), intent(in), optional :: panel_title
        type(plot_data_2d) :: estimated_curve
        type(plot_data_2d) :: exact_curve
        class(legend), pointer :: key

        ! The multiplot shares subplot axes, so keep them alive through draw.
        associate (axis_plot => axes(row,col))
            call axis_plot%initialize()
            call axis_plot%set_x_axis_title("Frequency [Hz]")
            call axis_plot%set_y_axis_title(label)
            if (present(panel_title)) call axis_plot%set_title(trim(panel_title))
            call estimated_curve%define_data(frequency, estimated_values)
            call estimated_curve%set_name("Estimated")
            call exact_curve%define_data(frequency, exact_values)
            call exact_curve%set_name("Analytical")
            call exact_curve%set_line_style(LINE_DASHED)
            call axis_plot%push(estimated_curve)
            call axis_plot%push(exact_curve)
            if (row == 1 .and. col == 1) then
                key => axis_plot%get_legend()
                call key%set_is_visible(.true.)
            end if
            call plt%set(row, col, axis_plot)
        end associate
    end subroutine
end program
