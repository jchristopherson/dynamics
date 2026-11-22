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
module dynamics_controls_tuning
    !! Defines an integral state feedback controller along with tools to
    !! evaluate its integral squared error (ISE) in closed form and to optimize
    !! its gains.  Gradients of the ISE with respect to the gains are computed
    !! with automatic differentiation.
    use iso_fortran_env, only : real64
    use ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use autodiff, only : dual, scalar_function, operator(+), operator(-), &
        operator(*), operator(/), abs
    use nonlin, only : bfgs, fcnnvar_helper
    use linalg, only : eigen
    use dynamics_error_handling, only : DYN_INVALID_INPUT_ERROR
    use dynamics_controls, only : state_space
    use dynamics_discrete_controls, only : discrete_state_space
    use dynamics_helper, only : solve_linear
    implicit none
    private
    public :: integral_state_feedback_controller
    public :: controller_optimizer

    real(real64), parameter :: UNSTABLE_PENALTY = 1.0d12
        !! The cost assigned to candidate gains that produce an unstable
        !! closed-loop system.

! ------------------------------------------------------------------------------
    type plant_data
        !! A real-valued copy of a plant used internally by the tuning routines.
        !! A zero sample time denotes a continuous-time plant.
        real(real64), allocatable, dimension(:,:) :: A
        real(real64), allocatable, dimension(:,:) :: B
        real(real64), allocatable, dimension(:,:) :: C
        real(real64), allocatable, dimension(:,:) :: D
        real(real64) :: sample_time = 0.0d0
        real(real64) :: control_weight = 0.0d0
    end type

! ------------------------------------------------------------------------------
    type integral_state_feedback_controller
        !! Defines a state feedback controller with integral action of the
        !! form
        !! $$ u = -Kx - K_i \xi, \qquad \dot{\xi} = y - r $$
        !! where \(K\) is the M-by-N state gain, \(K_i\) is the M-by-P integral
        !! gain, \(r\) is the reference, and \(\xi\) is the integral of the
        !! tracking error.
        real(real64), allocatable, dimension(:,:) :: K
            !! The M-by-N state feedback gain matrix.
        real(real64), allocatable, dimension(:,:) :: Ki
            !! The M-by-P integral gain matrix.
    contains
        generic, public :: is_stable => is_stable_continuous, is_stable_discrete
        generic, public :: integral_squared_error => ise_continuous, ise_discrete
        procedure, private :: is_stable_continuous
        procedure, private :: is_stable_discrete
        procedure, private :: ise_continuous
        procedure, private :: ise_discrete
    end type

    interface integral_state_feedback_controller
        module procedure :: isfc_init
    end interface

! ------------------------------------------------------------------------------
    type controller_optimizer
        !! Defines the settings used to optimize the gains of an
        !! [[integral_state_feedback_controller]] by minimizing the integral
        !! squared error of the closed-loop response to a unit step in the
        !! reference, plus an optional penalty on the control effort.
        real(real64) :: tolerance = 1.0d-10
            !! The convergence tolerance of the optimizer.
        integer :: max_evaluations = 1000
            !! The maximum number of function evaluations.
        real(real64) :: control_weight = 0.0d0
            !! The non-negative weight applied to the control effort, measured
            !! as the integral of the squared rate of change of the control
            !! input (the first difference for discrete systems).  A value of
            !! zero minimizes the integral squared error alone.
    contains
        generic, public :: optimize => optimize_continuous, optimize_discrete
        procedure, private :: optimize_continuous
        procedure, private :: optimize_discrete
    end type

contains
! ******************************************************************************
! CONTROLLER PROCEDURES
! ------------------------------------------------------------------------------
function isfc_init(K, Ki) result(rst)
    !! Initializes an integral state feedback controller.
    real(real64), intent(in), dimension(:,:) :: K
        !! The M-by-N state feedback gain matrix.
    real(real64), intent(in), dimension(:,:) :: Ki
        !! The M-by-P integral gain matrix.
    type(integral_state_feedback_controller) :: rst
        !! The controller.

    if (size(K, 1) /= size(Ki, 1)) error stop DYN_INVALID_INPUT_ERROR
    rst%K = K
    rst%Ki = Ki
end function

! ------------------------------------------------------------------------------
function is_stable_continuous(this, plant) result(rst)
    !! Determines if the closed-loop system is stable for a continuous-time
    !! plant.
    class(integral_state_feedback_controller), intent(in) :: this
        !! The [[integral_state_feedback_controller]] object.
    type(state_space), intent(in) :: plant
        !! The continuous-time plant.
    logical :: rst
        !! Returns true if the closed-loop system is stable; else, false.
    rst = plant_is_stable(plant_from_continuous(plant), this%K, this%Ki)
end function

! ------------------------------------------------------------------------------
function is_stable_discrete(this, plant) result(rst)
    !! Determines if the closed-loop system is stable for a discrete-time
    !! plant.
    class(integral_state_feedback_controller), intent(in) :: this
        !! The [[integral_state_feedback_controller]] object.
    type(discrete_state_space), intent(in) :: plant
        !! The discrete-time plant.
    logical :: rst
        !! Returns true if the closed-loop system is stable; else, false.
    rst = plant_is_stable(plant_from_discrete(plant), this%K, this%Ki)
end function

! ------------------------------------------------------------------------------
function ise_continuous(this, plant) result(rst)
    !! Computes the integral squared error of the closed-loop response to a
    !! unit step in the reference for a continuous-time plant.  An unstable
    !! closed loop returns positive infinity.
    class(integral_state_feedback_controller), intent(in) :: this
        !! The [[integral_state_feedback_controller]] object.
    type(state_space), intent(in) :: plant
        !! The continuous-time plant.
    real(real64) :: rst
        !! The integral squared error.
    rst = plant_ise(plant_from_continuous(plant), this%K, this%Ki)
end function

! ------------------------------------------------------------------------------
function ise_discrete(this, plant) result(rst)
    !! Computes the integral squared error of the closed-loop response to a
    !! unit step in the reference for a discrete-time plant.  The sum of
    !! squared errors is scaled by the sample time to approximate the
    !! continuous-time integral.  An unstable closed loop returns positive
    !! infinity.
    class(integral_state_feedback_controller), intent(in) :: this
        !! The [[integral_state_feedback_controller]] object.
    type(discrete_state_space), intent(in) :: plant
        !! The discrete-time plant.
    real(real64) :: rst
        !! The integral squared error.
    rst = plant_ise(plant_from_discrete(plant), this%K, this%Ki)
end function

! ******************************************************************************
! OPTIMIZER PROCEDURES
! ------------------------------------------------------------------------------
subroutine optimize_continuous(this, plant, controller, fout)
    !! Optimizes the gains of an integral state feedback controller for a
    !! continuous-time plant.  The controller must be stable on input.
    class(controller_optimizer), intent(in) :: this
        !! The [[controller_optimizer]] object.
    type(state_space), intent(in) :: plant
        !! The continuous-time plant.
    type(integral_state_feedback_controller), intent(inout) :: controller
        !! On input, the initial gains.  On output, the optimized gains.
    real(real64), intent(out), optional :: fout
        !! An optional output containing the optimized cost: the integral
        !! squared error plus the weighted control rate.
    call run_optimizer(this, plant_from_continuous(plant), controller, fout)
end subroutine

! ------------------------------------------------------------------------------
subroutine optimize_discrete(this, plant, controller, fout)
    !! Optimizes the gains of an integral state feedback controller for a
    !! discrete-time plant.  The controller must be stable on input.
    class(controller_optimizer), intent(in) :: this
        !! The [[controller_optimizer]] object.
    type(discrete_state_space), intent(in) :: plant
        !! The discrete-time plant.
    type(integral_state_feedback_controller), intent(inout) :: controller
        !! On input, the initial gains.  On output, the optimized gains.
    real(real64), intent(out), optional :: fout
        !! An optional output containing the optimized cost: the integral
        !! squared error plus the weighted control rate.
    call run_optimizer(this, plant_from_discrete(plant), controller, fout)
end subroutine

! ------------------------------------------------------------------------------
subroutine run_optimizer(this, plant, controller, fout)
    !! Packs the gains into a vector, and minimizes the integral squared
    !! error using BFGS with an AD gradient.
    class(controller_optimizer), intent(in) :: this
        !! The [[controller_optimizer]] object.
    type(plant_data), intent(in) :: plant
        !! The plant.
    type(integral_state_feedback_controller), intent(inout) :: controller
        !! The controller.
    real(real64), intent(out), optional :: fout
        !! An optional output containing the integral squared error at the
        !! optimized gains.

    ! Local Variables
    type(bfgs) :: solver
    type(fcnnvar_helper) :: obj
    type(plant_data) :: pd
    procedure(scalar_function), pointer :: fcn
    real(real64), allocatable, dimension(:) :: x
    integer :: m, n, p, nk

    if (this%control_weight < 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    call check_gain_dimensions(plant, controller%K, controller%Ki)
    if (.not.plant_is_stable(plant, controller%K, controller%Ki)) then
        error stop DYN_INVALID_INPUT_ERROR
    end if
    pd = plant
    pd%control_weight = this%control_weight

    m = size(controller%K, 1)
    n = size(controller%K, 2)
    p = size(controller%Ki, 2)
    nk = m * n
    allocate(x(nk + m * p))
    x(1:nk) = reshape(controller%K, [nk])
    x(nk + 1:) = reshape(controller%Ki, [m * p])

    fcn => isf_objective
    call obj%set_auto_fcn(fcn, size(x))
    call solver%set_max_fcn_evals(this%max_evaluations)
    call solver%set_tolerance(this%tolerance)
    call solver%set_use_line_search(.true.)
    call solver%solve(obj, x, fout, args = pd)

    controller%K = reshape(x(1:nk), [m, n])
    controller%Ki = reshape(x(nk + 1:), [m, p])
end subroutine

! ------------------------------------------------------------------------------
function isf_objective(x, args) result(f)
    !! The optimization objective.  The gains are unpacked from x, and the
    !! ISE is evaluated in dual arithmetic so that its derivatives are
    !! propagated by AD.  Unstable candidates receive a constant penalty with
    !! zero derivative.
    type(dual), intent(in), dimension(:) :: x
        !! The packed gains [vec(K); vec(Ki)].
    class(*), intent(in), optional :: args
        !! A [[plant_data]] object describing the plant.
    type(dual) :: f
        !! The cost.

    ! Local Variables
    integer :: m, n, p, nk
    type(dual), allocatable, dimension(:,:) :: K, Ki, Acl, Ce
    real(real64), allocatable, dimension(:,:) :: Bcl

    select type (args)
    type is (plant_data)
        m = size(args%B, 2)
        n = size(args%A, 1)
        p = size(args%C, 1)
        nk = m * n
        K = reshape(x(1:nk), [m, n])
        Ki = reshape(x(nk + 1:nk + m * p), [m, p])
        call build_closed_loop(args, K, Ki, Acl, Bcl, Ce)
        if (is_stable_loop(dual_values(Acl), args%sample_time)) then
            f = ise_from_loop(args, Acl, Bcl, Ce, gain_matrix(K, Ki))
        else
            f = dual(UNSTABLE_PENALTY, 0.0d0)
        end if
    class default
        error stop DYN_INVALID_INPUT_ERROR
    end select
end function

! ******************************************************************************
! PLANT CONVERSION
! ------------------------------------------------------------------------------
function plant_from_continuous(sys) result(rst)
    !! Creates a [[plant_data]] object from a continuous-time state-space model.
    type(state_space), intent(in) :: sys
        !! The continuous-time model.
    type(plant_data) :: rst
        !! The plant data.
    rst%A = sys%A
    rst%B = sys%B
    rst%C = sys%C
    rst%D = sys%D
    rst%sample_time = 0.0d0
end function

! ------------------------------------------------------------------------------
function plant_from_discrete(sys) result(rst)
    !! Creates a [[plant_data]] object from a discrete-time state-space model.
    type(discrete_state_space), intent(in) :: sys
        !! The discrete-time model.
    type(plant_data) :: rst
        !! The plant data.
    rst%A = sys%A
    rst%B = sys%B
    rst%C = sys%C
    rst%D = sys%D
    rst%sample_time = sys%sample_time
end function

! ------------------------------------------------------------------------------
subroutine check_gain_dimensions(plant, K, Ki)
    !! Verifies the gain matrices are consistent with the plant.
    type(plant_data), intent(in) :: plant
        !! The plant.
    real(real64), intent(in), dimension(:,:) :: K
        !! The state feedback gain.
    real(real64), intent(in), dimension(:,:) :: Ki
        !! The integral gain.

    ! Local Variables
    integer :: m, n, p

    m = size(plant%B, 2)
    n = size(plant%A, 1)
    p = size(plant%C, 1)
    if (size(K, 1) /= m .or. size(K, 2) /= n) error stop DYN_INVALID_INPUT_ERROR
    if (size(Ki, 1) /= m .or. size(Ki, 2) /= p) error stop DYN_INVALID_INPUT_ERROR
end subroutine

! ******************************************************************************
! CLOSED-LOOP ANALYSIS
! ------------------------------------------------------------------------------
function plant_is_stable(plant, K, Ki) result(rst)
    !! Determines if the closed-loop system formed by the plant and the gains
    !! is stable.
    type(plant_data), intent(in) :: plant
        !! The plant.
    real(real64), intent(in), dimension(:,:) :: K
        !! The state feedback gain.
    real(real64), intent(in), dimension(:,:) :: Ki
        !! The integral gain.
    logical :: rst
        !! Returns true if the closed loop is stable; else, false.

    ! Local Variables
    type(dual), allocatable, dimension(:,:) :: Kd, Kid, Acl, Ce
    real(real64), allocatable, dimension(:,:) :: Bcl

    call check_gain_dimensions(plant, K, Ki)
    Kd = to_dual_matrix(K)
    Kid = to_dual_matrix(Ki)
    call build_closed_loop(plant, Kd, Kid, Acl, Bcl, Ce)
    rst = is_stable_loop(dual_values(Acl), plant%sample_time)
end function

! ------------------------------------------------------------------------------
function plant_ise(plant, K, Ki) result(rst)
    !! Computes the integral squared error for the given plant and gains.
    type(plant_data), intent(in) :: plant
        !! The plant.
    real(real64), intent(in), dimension(:,:) :: K
        !! The state feedback gain.
    real(real64), intent(in), dimension(:,:) :: Ki
        !! The integral gain.
    real(real64) :: rst
        !! The integral squared error, or positive infinity if unstable.

    ! Local Variables
    type(dual), allocatable, dimension(:,:) :: Kd, Kid, Acl, Ce
    real(real64), allocatable, dimension(:,:) :: Bcl
    type(dual) :: cost

    call check_gain_dimensions(plant, K, Ki)
    Kd = to_dual_matrix(K)
    Kid = to_dual_matrix(Ki)
    call build_closed_loop(plant, Kd, Kid, Acl, Bcl, Ce)
    if (.not.is_stable_loop(dual_values(Acl), plant%sample_time)) then
        rst = ieee_value(rst, ieee_positive_inf)
        return
    end if
    cost = ise_from_loop(plant, Acl, Bcl, Ce, gain_matrix(Kd, Kid))
    rst = cost%x
end function

! ------------------------------------------------------------------------------
subroutine build_closed_loop(plant, K, Ki, Acl, Bcl, Ce)
    !! Builds the closed-loop system for the augmented state z = [x; xi], where
    !! xi is the integral of the tracking error, so xi' = y - r.  For a discrete-time plant,
    !! the integrator is updated with the sample time, Ts.  The error is
    !! e = y - r = Ce z.
    type(plant_data), intent(in) :: plant
        !! The plant.
    type(dual), intent(in), dimension(:,:) :: K
        !! The state feedback gain.
    type(dual), intent(in), dimension(:,:) :: Ki
        !! The integral gain.
    type(dual), intent(out), allocatable, dimension(:,:) :: Acl
        !! The closed-loop dynamics matrix.
    real(real64), intent(out), allocatable, dimension(:,:) :: Bcl
        !! The closed-loop reference input matrix.
    type(dual), intent(out), allocatable, dimension(:,:) :: Ce
        !! The error output matrix.

    ! Local Variables
    integer :: n, m, p, nz, i, j, l
    real(real64) :: h, d
    type(dual) :: s

    n = size(plant%A, 1)
    m = size(plant%B, 2)
    p = size(plant%C, 1)
    nz = n + p
    if (plant%sample_time > 0.0d0) then
        h = plant%sample_time
        d = 1.0d0
    else
        h = 1.0d0
        d = 0.0d0
    end if

    allocate(Acl(nz, nz), Bcl(nz, p), Ce(p, nz))
    Bcl = 0.0d0
    do i = 1, p
        Bcl(n + i, i) = -h
    end do

    ! Top-left: A - B*K
    do j = 1, n
        do i = 1, n
            s = dual(plant%A(i, j), 0.0d0)
            do l = 1, m
                s = s - plant%B(i, l) * K(l, j)
            end do
            Acl(i, j) = s
        end do
    end do

    ! Top-right: -B*Ki
    do j = 1, p
        do i = 1, n
            s = dual(0.0d0, 0.0d0)
            do l = 1, m
                s = s - plant%B(i, l) * Ki(l, j)
            end do
            Acl(i, n + j) = s
        end do
    end do

    ! Bottom-left: h*(C - D*K)
    do j = 1, n
        do i = 1, p
            s = dual(h * plant%C(i, j), 0.0d0)
            do l = 1, m
                s = s - h * plant%D(i, l) * K(l, j)
            end do
            Acl(n + i, j) = s
        end do
    end do

    ! Bottom-right: d*I - h*D*Ki
    do j = 1, p
        do i = 1, p
            if (i == j) then
                s = dual(d, 0.0d0)
            else
                s = dual(0.0d0, 0.0d0)
            end if
            do l = 1, m
                s = s - h * plant%D(i, l) * Ki(l, j)
            end do
            Acl(n + i, n + j) = s
        end do
    end do

    ! Error output: e = (C - D*K)x - D*Ki*xi - r
    do j = 1, n
        do i = 1, p
            s = dual(plant%C(i, j), 0.0d0)
            do l = 1, m
                s = s - plant%D(i, l) * K(l, j)
            end do
            Ce(i, j) = s
        end do
    end do
    do j = 1, p
        do i = 1, p
            s = dual(0.0d0, 0.0d0)
            do l = 1, m
                s = s - plant%D(i, l) * Ki(l, j)
            end do
            Ce(i, n + j) = s
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
function ise_from_loop(plant, Acl, Bcl, Ce, G) result(cost)
    !! Computes the ISE for a unit step in every reference channel, plus the
    !! weighted control rate when plant%control_weight is positive.  The
    !! initial state is zero, so z(0) - z_ss = -z_ss.  The gramian P satisfies
    !! a Lyapunov equation (continuous) or a Stein equation (discrete), and the
    !! cost is z0' P z0.  The control input is u - u_ss = -G (z - z_ss), where
    !! G = [K, Ki] and z = [x; xi], so its rate uses the deviation dynamics
    !! z - z_ss evolving under Acl.
    type(plant_data), intent(in) :: plant
        !! The plant.
    type(dual), intent(in), dimension(:,:) :: Acl
        !! The closed-loop dynamics matrix.
    real(real64), intent(in), dimension(:,:) :: Bcl
        !! The closed-loop reference input matrix.
    type(dual), intent(in), dimension(:,:) :: Ce
        !! The error output matrix.
    type(dual), intent(in), dimension(:,:) :: G
        !! The M-by-(N+P) gain matrix [K, Ki] acting on the augmented state.
    type(dual) :: cost
        !! The integral squared error plus the control effort penalty.

    ! Local Variables
    integer :: nz, np, i, j, k, a, b, r, c
    logical :: discrete
    real(real64) :: scale
    type(dual) :: s, dm
    type(dual), allocatable, dimension(:,:) :: Ms, Q, M, H
    type(dual), allocatable, dimension(:) :: z, rhs

    nz = size(Acl, 1)
    np = size(Ce, 1)
    discrete = plant%sample_time > 0.0d0

    ! Steady-state response, z_ss, to a unit step in every reference
    allocate(Ms(nz, nz), z(nz))
    do j = 1, nz
        do i = 1, nz
            Ms(i, j) = -Acl(i, j)
        end do
        if (discrete) Ms(j, j) = Ms(j, j) + 1.0d0
    end do
    do i = 1, nz
        z(i) = dual(sum(Bcl(i, :)), 0.0d0)
    end do
    call solve_linear(Ms, z)
    do i = 1, nz
        z(i) = -z(i)
    end do

    ! Gramian input weight, scaled by the sample time for discrete systems
    if (discrete) then
        scale = plant%sample_time
    else
        scale = 1.0d0
    end if
    allocate(Q(nz, nz))
    do j = 1, nz
        do i = 1, nz
            s = dual(0.0d0, 0.0d0)
            do k = 1, np
                s = s + Ce(k, i) * Ce(k, j)
            end do
            Q(i, j) = scale * s
        end do
    end do
    if (plant%control_weight > 0.0d0) then
        ! The control rate, du/dt = -G Acl (z - z_ss), is penalized.  For
        ! discrete systems the first difference, -G (Acl - I) (z - z_ss), is
        ! used, and the sum is divided by the sample time to approximate the
        ! integral.  Penalizing the rate, rather than the magnitude, does not
        ! reward faster loops.
        allocate(H(size(G, 1), nz))
        do j = 1, nz
            do i = 1, size(G, 1)
                s = dual(0.0d0, 0.0d0)
                do k = 1, nz
                    dm = Acl(k, j)
                    if (discrete .and. k == j) dm = dm - 1.0d0
                    s = s + G(i, k) * dm
                end do
                H(i, j) = s
            end do
        end do
        do j = 1, nz
            do i = 1, nz
                s = dual(0.0d0, 0.0d0)
                do k = 1, size(G, 1)
                    s = s + H(k, i) * H(k, j)
                end do
                Q(i, j) = Q(i, j) + plant%control_weight * s / scale
            end do
        end do
    end if

    ! Solve for the gramian, vec(P), using the Kronecker form
    allocate(M(nz * nz, nz * nz), rhs(nz * nz))
    M = dual(0.0d0, 0.0d0)
    do j = 1, nz
        do i = 1, nz
            r = (j - 1) * nz + i
            if (discrete) then
                M(r, r) = M(r, r) + 1.0d0
                do b = 1, nz
                    do a = 1, nz
                        c = (b - 1) * nz + a
                        M(r, c) = M(r, c) - Acl(a, i) * Acl(b, j)
                    end do
                end do
                rhs(r) = Q(i, j)
            else
                do k = 1, nz
                    c = (j - 1) * nz + k
                    M(r, c) = M(r, c) + Acl(k, i)
                    c = (k - 1) * nz + i
                    M(r, c) = M(r, c) + Acl(k, j)
                end do
                rhs(r) = -Q(i, j)
            end if
        end do
    end do
    call solve_linear(M, rhs)

    cost = dual(0.0d0, 0.0d0)
    do j = 1, nz
        do i = 1, nz
            cost = cost + z(i) * rhs((j - 1) * nz + i) * z(j)
        end do
    end do
end function



! ------------------------------------------------------------------------------
pure function is_stable_loop(Acl, sample_time) result(rst)
    !! Determines if a closed-loop dynamics matrix is stable.  Continuous
    !! systems require all eigenvalues to lie in the open left half plane.
    !! Discrete systems require all eigenvalues to lie inside the unit circle.
    real(real64), intent(in), dimension(:,:) :: Acl
        !! The closed-loop dynamics matrix.
    real(real64), intent(in) :: sample_time
        !! The sample time, or zero for continuous-time systems.
    logical :: rst
        !! Returns true if stable; else, false.

    ! Local Variables
    complex(real64), allocatable, dimension(:) :: lambda

    allocate(lambda(size(Acl, 1)))
    call eigen(Acl, lambda)
    if (sample_time > 0.0d0) then
        rst = all(abs(lambda) < 1.0d0)
    else
        rst = all(real(lambda) < 0.0d0)
    end if
end function

! ------------------------------------------------------------------------------
function gain_matrix(K, Ki) result(rst)
    !! Concatenates the state feedback and integral gains into [K, Ki], which
    !! maps the augmented state z = [x; xi] to the control input.
    type(dual), intent(in), dimension(:,:) :: K
        !! The M-by-N state feedback gain.
    type(dual), intent(in), dimension(:,:) :: Ki
        !! The M-by-P integral gain.
    type(dual), allocatable, dimension(:,:) :: rst
        !! The M-by-(N+P) gain matrix.

    integer :: n

    n = size(K, 2)
    allocate(rst(size(K, 1), n + size(Ki, 2)))
    rst(:, 1:n) = K
    rst(:, n + 1:) = Ki
end function

! ------------------------------------------------------------------------------
pure function to_dual_matrix(a) result(rst)
    !! Converts a real matrix to a dual matrix with zero derivatives.
    real(real64), intent(in), dimension(:,:) :: a
        !! The real matrix.
    type(dual), allocatable, dimension(:,:) :: rst
        !! The dual matrix.

    ! Local Variables
    integer :: i, j

    allocate(rst(size(a, 1), size(a, 2)))
    do j = 1, size(a, 2)
        do i = 1, size(a, 1)
            rst(i, j) = dual(a(i, j), 0.0d0)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function dual_values(a) result(rst)
    !! Extracts the values of a dual matrix.
    type(dual), intent(in), dimension(:,:) :: a
        !! The dual matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The real values.
    allocate(rst(size(a, 1), size(a, 2)))
    rst = a%x
end function

! ------------------------------------------------------------------------------
end module