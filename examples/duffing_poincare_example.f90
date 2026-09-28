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
module duffing_equation
    use iso_fortran_env, only : real64
    implicit none

    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 5.0d0
    real(real64), parameter :: delta = 2.0d-2
    real(real64), parameter :: gamma = 8.0d0
    real(real64), parameter :: omega = 5.0d-1

contains

subroutine duffing_eom(t, state, derivative, args)
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: state
    real(real64), intent(out), dimension(:) :: derivative
    class(*), intent(inout), optional :: args

    ! State is [displacement, velocity, forcing phase]. The phase state makes
    ! the periodically forced equation autonomous for section construction.
    derivative(1) = state(2)
    derivative(2) = gamma * cos(state(3)) - delta * state(2) - &
        alpha * state(1) - beta * state(1)**3
    derivative(3) = omega
end subroutine

subroutine duffing_section_coordinates(t, state, coordinates_out)
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: state
    real(real64), intent(out), dimension(3) :: coordinates_out

    coordinates_out = [state(1), state(2), sin(state(3))]
end subroutine

end module

program duffing_poincare_example
    use iso_fortran_env, only : int32, real64
    use duffing_equation
    use diffeq
    use dynamics_maps, only : poincare_map, POINCARE_ONE_SIDED_FROM_BACK
    use fplot_core
    implicit none

    integer(int32), parameter :: transient_cycles = 50
    integer(int32), parameter :: response_cycles = 50000
    integer(int32), parameter :: points_per_cycle = 50
    integer(int32), parameter :: response_points = &
        response_cycles * points_per_cycle + 1
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: period = 2.0d0 * pi / omega
    real(real64), parameter :: transient_end = transient_cycles * period
    real(real64), parameter :: response_end = &
        (transient_cycles + response_cycles) * period

    real(real64), allocatable, dimension(:,:) :: transient_solution, section
    type(ode_container) :: model
    type(runge_kutta_45) :: integrator
    type(plot_2d) :: section_plot
    type(plot_data_2d) :: section_data

    ! Solve the transient using only its endpoint as the response initial state.
    model%fcn => duffing_eom
    call integrator%solve(model, [0.0d0, 0.5d0 * transient_end, &
        transient_end], [0.0d0, 0.0d0, 0.0d0])
    transient_solution = integrator%get_solution()

    ! In the autonomous embedding, sin(phase)=0 is crossed from back to front
    ! once per forcing period at phase 0 modulo 2*pi. Retain only the section.
    section = poincare_map(model, [transient_end, response_end], &
        transient_solution(3,2:4), response_points, &
        side = POINCARE_ONE_SIDED_FROM_BACK, solver = integrator, &
        chunk_size = 20 * points_per_cycle, &
        coordinates = duffing_section_coordinates)

    ! Plot the Poincare points in displacement-velocity coordinates.
    call section_plot%initialize()
    call section_plot%set_x_axis_title("x")
    call section_plot%set_y_axis_title("dx/dt")
    call section_data%define_data(section(:,1), section(:,2))
    call section_data%set_draw_line(.false.)
    call section_data%set_draw_markers(.true.)
    call section_data%set_marker_scaling(0.5)
    call section_data%set_line_color(CLR_BLACK)
    call section_plot%push(section_data)
    call section_plot%draw()
end program
