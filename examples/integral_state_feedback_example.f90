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
module step_model
    use iso_fortran_env
    implicit none
contains
    subroutine step_command(t, r, args)
        !! Returns a unit step in the reference, applied at t = 0.
        real(real64), intent(in) :: t
            !! The current simulation time
        real(real64), intent(out) :: r(:)
            !! The reference signal.
        class(*), intent(in), optional :: args
            !! Unused.
        r = 1.0d0
    end subroutine
end module

program integral_state_feedback_example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use step_model
    implicit none

    ! Solution parameters
    real(real64), parameter :: tmax = 10.0d0

    ! Local variables
    integer :: i, n
    real(real64) :: t(2), ic(3), ise_cf, ise_num, fout
    real(real64), allocatable, dimension(:) :: r
    real(real64), allocatable, dimension(:,:) :: A, B, C, D, K, Ki
    real(real64), allocatable, dimension(:,:) :: sol0, sol1, sol2
    type(state_space) :: plant, mdl0, mdl1, mdl2
    type(integral_state_feedback_controller) :: ctrl, ctrl_pen
    type(controller_optimizer) :: opt

    ! Plot variables
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis, x2, y2
    type(multiplot) :: fig
    class(legend), pointer :: lgnd

    ! Define a continuous-time plant with a feedthrough term.  The state is
    ! [x1; x2] and the input is u.
    allocate(A(2,2), B(2,1), C(1,2), D(1,1))
    A = reshape([0.0d0, -2.0d0, 1.0d0, -3.0d0], [2, 2])
    B = reshape([0.0d0, 1.0d0], [2, 1])
    C = reshape([1.0d0, 0.0d0], [1, 2])
    D = reshape([0.2d0], [1, 1])
    plant = state_space(A, B, C, D)
    n = size(A, 1)

    ! Define a stable starting controller
    K = reshape([2.0d0, 1.0d0], [1, 2])
    Ki = reshape([1.0d0], [1, 1])
    ctrl = integral_state_feedback_controller(K, Ki)

    ! Closed-form ISE and stability of the starting controller
    print "(A, L1)", "Stable (initial gains): ", ctrl%is_stable(plant)
    ise_cf = ctrl%integral_squared_error(plant)
    print "(A, F0.6)", "ISE (initial gains): ", ise_cf

    ! Optimize the gains.  The gradient of the ISE is computed with automatic
    ! differentiation.
    opt%tolerance = 1.0d-8
    opt%max_evaluations = 1000
    call opt%optimize(plant, ctrl, fout)

    print "(A, L1)", "Stable (optimized gains): ", ctrl%is_stable(plant)
    print "(A, F0.6)", "ISE (optimized gains): ", ctrl%integral_squared_error(plant)
    print "(A, F0.6)", "ISE reported by optimizer: ", fout
    print "(A)", "Optimized K:"
    print "(A, 2F0.6)", achar(9), ctrl%K
    print "(A)", "Optimized Ki:"
    print "(A, F0.6)", achar(9), ctrl%Ki(1,1)

    ! Simulate the closed-loop response to a unit step with the optimized
    ! gains and compare the numerical ISE to the closed-form value.
    t = [0.0d0, tmax]
    ic = 0.0d0
    mdl1 = ctrl%apply_to_state_space_plant(plant)
    sol1 = lti_solve(mdl1, step_command, t, ic)
    ise_num = integral_error(sol1)
    print "(A, F0.6)", "Numerical ISE (optimized gains): ", ise_num
    print "(A, I0)", "Solution points (optimized gains): ", size(sol1, 1)

    ! Simulate the response with the starting gains for comparison
    ctrl = integral_state_feedback_controller(K, Ki)
    mdl0 = ctrl%apply_to_state_space_plant(plant)
    sol0 = lti_solve(mdl0, step_command, t, ic)
    print "(A, F0.6)", "Numerical ISE (initial gains): ", integral_error(sol0)

    ! Optimize again from the starting gains with a control effort weight.  The
    ! optimizer then minimizes the ISE plus the weighted control effort, which
    ! limits the size of the gains.
    ctrl_pen = integral_state_feedback_controller(K, Ki)
    opt%control_weight = 1.0d0
    call opt%optimize(plant, ctrl_pen, fout)

    print "(A, L1)", "Stable (control weight, optimized gains): ", ctrl_pen%is_stable(plant)
    print "(A, F0.6)", "ISE (control weight, optimized gains): ", ctrl_pen%integral_squared_error(plant)
    print "(A, F0.6)", "Cost reported by optimizer (control weight): ", fout
    print "(A)", "Control weight optimized K:"
    print "(A, 2F0.6)", achar(9), ctrl_pen%K
    print "(A)", "Control weight optimized Ki:"
    print "(A, F0.6)", achar(9), ctrl_pen%Ki(1,1)
    mdl2 = ctrl_pen%apply_to_state_space_plant(plant)
    sol2 = lti_solve(mdl2, step_command, t, ic)
    print "(A, F0.6)", "Numerical ISE (control weight, optimized gains): ", integral_error(sol2)

    ! Plot the reference and output (top) and control effort (bottom)
    allocate(r(size(sol1, 1)))
    r = 1.0d0
    call fig%initialize(2, 1, width = 1000, height = 600)
    call plt1%initialize()
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    lgnd => plt1%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call xAxis%set_title("t [s]")
    call yAxis%set_title("y(t)")

    call pd%define_data(sol0(:,1), sol0(:,3))
    call pd%set_line_width(2.0)
    call pd%set_name("Initial gains")
    call plt1%push(pd)

    call pd%clear_all()
    call pd%define_data(sol1(:,1), sol1(:,3))
    call pd%set_line_width(2.0)
    call pd%set_name("Optimized gains")
    call plt1%push(pd)

    call pd%clear_all()
    call pd%define_data(sol2(:,1), sol2(:,3))
    call pd%set_line_width(2.0)
    call pd%set_name("Control weight")
    call plt1%push(pd)

    call pd%clear_all()
    call pd%define_data(sol1(:,1), r)
    call pd%set_name("Reference")
    call plt1%push(pd)

    call plt2%initialize()
    x2 => plt2%get_x_axis()
    y2 => plt2%get_y_axis()
    call x2%set_title("t [s]")
    call y2%set_title("u(t)")

    call pd%clear_all()
    call pd%define_data(sol0(:,1), sol0(:,2))
    call pd%set_line_width(2.0)
    call pd%set_name("Initial gains")
    call plt2%push(pd)

    call pd%clear_all()
    call pd%define_data(sol1(:,1), sol1(:,2))
    call pd%set_line_width(2.0)
    call pd%set_name("Optimized gains")
    call plt2%push(pd)

    call pd%clear_all()
    call pd%define_data(sol2(:,1), sol2(:,2))
    call pd%set_line_width(2.0)
    call pd%set_name("Control weight")
    call plt2%push(pd)

    call fig%set(1, 1, plt1)
    call fig%set(2, 1, plt2)
    call fig%draw()

contains

    function integral_error(sol) result(rst)
        !! Numerically integrates the squared tracking error (y - 1)^2 over
        !! the solution with the trapezoidal rule.
        real(real64), intent(in), dimension(:,:) :: sol
            !! The solution from simulate_closed_loop.
        real(real64) :: rst
            !! The integral squared error.

        integer :: j
        real(real64) :: e1, e2

        rst = 0.0d0
        do j = 1, size(sol, 1) - 1
            e1 = sol(j, 3) - 1.0d0
            e2 = sol(j+1, 3) - 1.0d0
            rst = rst + 0.5d0 * (e1 * e1 + e2 * e2) * (sol(j+1, 1) - sol(j, 1))
        end do
    end function
end program
