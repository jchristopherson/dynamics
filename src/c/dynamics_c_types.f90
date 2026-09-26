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
module dynamics_c_types
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use diffeq
    use spectrum, only : window
    use dynamics_error_handling
    use nonlin
    implicit none

    interface
        subroutine c_vecfcn(nvar, neqn, x, f) bind(C, name = "c_vecfcn")
            use iso_c_binding
            integer(c_int), intent(in), value :: nvar
            integer(c_int), intent(in), value :: neqn
            real(c_double), intent(in) :: x(nvar)
            real(c_double), intent(out) :: f(neqn)
        end subroutine

        subroutine c_modal_excite(n, freq, f) bind(C, name = "c_modal_excite")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: freq
            complex(c_double), intent(out) :: f(n)
        end subroutine

        subroutine c_harmonic_ode(n, freq, t, x, dxdt) &
            bind(C, name = "c_harmonic_ode")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: freq
            real(c_double), intent(in), value :: t
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out) :: dxdt(n)
        end subroutine

        pure function c_window_function(n, bin) result(rst) &
            bind(C, name = "c_window_function")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            integer(c_int), intent(in), value :: bin
            real(c_double) :: rst
        end function

        subroutine c_constraint_equations(n, neqn, nparam, xg, fg, xc, p, &
            fc) bind(C, name = "c_constraint_equations")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            integer(c_int), intent(in), value :: neqn
            integer(c_int), intent(in), value :: nparam
            real(c_double), intent(in) :: xg(n)
            real(c_double), intent(in) :: fg(n)
            real(c_double), intent(in) :: xc(neqn)
            real(c_double), intent(in) :: p(nparam)
            real(c_double), intent(out) :: fc(neqn)
        end subroutine

        subroutine c_ode_fit(neqn, nparam, mdl, t, x, F, dxdt) &
            bind(C, name = "c_ode_fit")
            use iso_c_binding
            integer(c_int), intent(in), value :: neqn
            integer(c_int), intent(in), value :: nparam
            real(c_double), intent(in) :: mdl(nparam)
            real(c_double), intent(in), value :: t
            real(c_double), intent(in) :: x(neqn)
            real(c_double), intent(in), value :: F
            real(c_double), intent(out) :: dxdt(neqn)
        end subroutine

        subroutine c_ss_excitation(n, t, u) bind(C, name = "c_ss_excitation")
            use iso_c_binding
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: t
            real(c_double), intent(out) :: u(n)
        end subroutine

    end interface

    type c_vecfcn_container
        procedure(c_vecfcn), pointer, nopass :: fcn
    end type

    type c_modal_excite_container
        procedure(c_modal_excite), pointer, nopass :: fcn
    end type

    type c_harmonic_ode_container
        procedure(c_harmonic_ode), pointer, nopass :: fcn
    end type

    type c_siso_fit_container
        procedure(c_ode_fit), pointer, nopass :: odefcn
        procedure(c_constraint_equations), pointer, nopass :: constraints
    end type

    type c_ss_excitation_container
        procedure(c_ss_excitation), pointer, nopass :: fcn
    end type

    type, bind(C) :: c_iteration_behavior
        logical(c_bool) :: converge_on_chng
        logical(c_bool) :: converge_on_fcn
        logical(c_bool) :: converge_on_zero_diff
        integer(c_int) :: fcn_count
        integer(c_int) :: gradient_count
        integer(c_int) :: iter_count
        integer(c_int) :: jacobian_count
    end type

    type, bind(C) :: c_frequency_sweep_controls
        integer(c_int) :: cycle_count
        integer(c_int) :: transient_cycles
        integer(c_int) :: points_per_cycle
        logical(c_bool) :: frequency_in_hz
    end type

    type, bind(C) :: c_iteration_controls
        real(c_double) :: change_in_solution_tolerance
        real(c_double) :: gradient_tolerance
        real(c_double) :: iteration_improvement_tolerance
        real(c_double) :: residual_tolerance
        integer(c_int) :: max_function_evaluations
        integer(c_int) :: max_iteration_between_updates
        integer(c_int) :: max_iteration_count
    end type

    type, bind(C) :: c_regression_statistics
        real(c_double) :: confidence_interval
        real(c_double) :: probability
        real(c_double) :: standard_error
        real(c_double) :: t_statistic
    end type
    
    type, bind(C) :: c_lm_solver_options
        real(c_double) :: damping_decrease_factor
        real(c_double) :: damping_increase_factor
        real(c_double) :: finite_difference_step_size
        integer(c_int) :: method
    end type

    type, extends(window) :: c_window
        procedure(c_window_function), pointer, nopass :: fcn
    contains
        procedure, public :: evaluate => cw_eval
    end type

    type, bind(C) :: c_dynamic_system_measurement
        integer(c_int) :: npts
        type(c_ptr) :: input
        type(c_ptr) :: output
        type(c_ptr) :: t
    end type

    integer(c_int), parameter :: DYN_RUNGE_KUTTA_23 = 10
    integer(c_int), parameter :: DYN_RUNGE_KUTTA_45 = 11
    integer(c_int), parameter :: DYN_RUNGE_KUTTA_853 = 12
    integer(c_int), parameter :: DYN_ROSENBROCK = 13
    integer(c_int), parameter :: DYN_BDF = 14
    integer(c_int), parameter :: DYN_ADAMS = 15
    integer(c_int), parameter :: DYN_KENNEDY_CARPENTER_4 = 16
    integer(c_int), parameter :: DYN_KENNEDY_CARPENTER_5 = 17
    integer(c_int), parameter :: DYN_TSITOURAS_5 = 18

    type, bind(C) :: c_quaternion
        real(c_double) :: w
        real(c_double) :: x
        real(c_double) :: y
        real(c_double) :: z
    end type

    type, bind(C) :: c_plane
        real(c_double) :: a
        real(c_double) :: b
        real(c_double) :: c
        real(c_double) :: d
    end type

    type, bind(C) :: c_line
        real(c_double) :: r0(3)
        real(c_double) :: v(3)
    end type

    type, bind(C) :: c_plucker_line
        real(c_double) :: u(3)
        real(c_double) :: m(3)
    end type

    type, bind(C) :: c_coordinate_system
        real(c_double) :: origin(3)
        real(c_double) :: i(3)
        real(c_double) :: j(3)
        real(c_double) :: k(3)
    end type

    type, bind(C) :: c_dh_parameter_set
        real(c_double) :: link_length
        real(c_double) :: link_twist
        real(c_double) :: link_offset
        real(c_double) :: joint_angle
    end type

    type, bind(C) :: c_dh_table
        integer(c_int) :: count
        type(c_ptr) :: parameters   ! pointer to an array of c_dh_parameter_set items
    end type

    type, bind(C) :: c_binary_link
        real(c_double) :: link_length
        real(c_double) :: link_twist
        real(c_double) :: link_offset
        real(c_double) :: joint_angle
        integer(c_int) :: joint_type
        real(c_double) :: mass
        real(c_double) :: cg(3)
        real(c_double) :: inertia(9)
    end type

    type, bind(C) :: c_serial_linkage
        integer(c_int) :: link_count
        type(c_ptr) :: links    ! pointer to an array of c_binary_link items
    end type

    type, bind(C) :: c_mechanism_link
        integer(c_int) :: frame_count
        type(c_ptr) :: frames   ! pointer to a 4-by-4-by-frame_count array
        real(c_double) :: mass
        real(c_double) :: cg(3)
        real(c_double) :: inertia(9)
    end type

    type, bind(C) :: c_joint
        integer(c_int) :: joint_type
        integer(c_int) :: parent_link
        integer(c_int) :: parent_frame
        integer(c_int) :: child_link
        integer(c_int) :: child_frame
        logical(c_bool) :: actuated
    end type

    type :: c_mechanism_container
        ! The object referenced by an opaque mechanism handle.
        class(kinematic_mechanism), allocatable :: item
    end type

    type, bind(C) :: c_polynomial
        integer(c_int) :: order
        type(c_ptr) :: coefficients ! ascending order
    end type

    type, bind(C) :: c_transfer_function
        type(c_polynomial) :: numerator
        type(c_polynomial) :: denominator
    end type

    type, bind(C) :: c_state_space_model
        integer(c_int) :: dimension
        integer(c_int) :: n_inputs
        integer(c_int) :: n_outputs
        type(c_ptr) :: A    ! dimension -by- dimension
        type(c_ptr) :: B    ! dimension -by- n_inputs
        type(c_ptr) :: C    ! n_outputs -by- dimension
        type(c_ptr) :: D    ! n_outputs -by- n_inputs
    end type

    type, bind(C) :: c_material
        real(c_double) :: density
        real(c_double) :: modulus
        real(c_double) :: poissons_ratio
    end type

    type, bind(C) :: c_node
        integer(c_int) :: index
        integer(c_int) :: dof
        real(c_double) :: x
        real(c_double) :: y
        real(c_double) :: z
    end type

    type, bind(C) :: c_beam_element_2d
        type(c_material) :: material
        real(c_double) :: area
        real(c_double) :: moment_of_inertia
        type(c_node) :: node_1
        type(c_node) :: node_2
    end type

    type, bind(C) :: c_beam_element_3d
        type(c_material) :: material
        real(c_double) :: area
        real(c_double) :: Ixx
        real(c_double) :: Iyy
        real(c_double) :: Izz
        real(c_double) :: Iyz
        type(c_node) :: node_1
        type(c_node) :: node_2
        real(c_double) :: orientation_point(3)
    end type

    interface assignment(=)
        module procedure :: convert_to_c_iteration_behavior
        module procedure :: convert_from_c_iteration_behavior
        module procedure :: convert_to_c_quaternion
        module procedure :: convert_from_c_quaternion
        module procedure :: convert_to_c_line
        module procedure :: convert_from_c_line
        module procedure :: convert_to_c_plane
        module procedure :: convert_from_c_plane
        module procedure :: convert_to_c_plucker_line
        module procedure :: convert_from_c_plucker_line
        module procedure :: convert_to_c_coordinate_system
        module procedure :: convert_from_c_coordinate_system
        module procedure :: convert_to_c_dh_parameter_set
        module procedure :: convert_from_c_dh_parameter_set
        module procedure :: convert_from_c_dh_table
        module procedure :: convert_to_c_binary_link
        module procedure :: convert_from_c_binary_link
        module procedure :: convert_from_c_serial_linkage
        module procedure :: convert_to_c_polynomial
        module procedure :: convert_from_c_polynomial
        module procedure :: convert_to_c_transfer_function
        module procedure :: convert_from_c_transfer_function
        module procedure :: convert_to_c_state_space
        module procedure :: convert_from_c_state_space
    end interface

    interface
        function c_alloc_dh_table(n, tbl) result(rst) &
            bind(C, name = "c_alloc_dh_table")
            use iso_c_binding, only : c_int
            import c_dh_table
            integer(c_int), intent(in), value :: n
            type(c_dh_table), intent(out) :: tbl
            integer(c_int) :: rst
        end function

        function c_alloc_serial_linkage(n, lnk) result(rst) &
            bind(C, name = "c_alloc_serial_linkage")
            use iso_c_binding, only : c_int
            import c_serial_linkage
            integer(c_int), intent(in), value :: n
            type(c_serial_linkage), intent(out) :: lnk
            integer(c_int) :: rst
        end function

        function c_alloc_polynomial(order, p) result(rst) &
            bind(C, name = "c_alloc_polynomial")
            use iso_c_binding, only : c_int
            import c_polynomial
            integer(c_int), intent(in), value :: order
            type(c_polynomial), intent(out) :: p
            integer(c_int) :: rst
        end function

        function c_alloc_transfer_function(numer_order, denom_order, tf) &
            result(rst) bind(C, name = "c_alloc_transfer_function")
            use iso_c_binding, only : c_int
            import c_transfer_function
            integer(c_int), intent(in), value :: numer_order, denom_order
            type(c_transfer_function), intent(out) :: tf
            integer(c_int) :: rst
        end function

        function c_alloc_state_space_model(dimension, n_inputs, n_outputs, &
            mdl) result(rst) bind(C, name = "c_alloc_state_space_model")
            use iso_c_binding, only : c_int
            import c_state_space_model
            integer(c_int), intent(in), value :: dimension
            integer(c_int), intent(in), value :: n_inputs
            integer(c_int), intent(in), value :: n_outputs
            type(c_state_space_model), intent(out) :: mdl
            integer(c_int) :: rst
        end function
    end interface

contains

pure subroutine convert_to_c_iteration_behavior(c, f)
    type(c_iteration_behavior), intent(out) :: c
    type(iteration_behavior), intent(in) :: f
    c%converge_on_chng = logical(f%converge_on_chng, c_bool)
    c%converge_on_fcn = logical(f%converge_on_fcn, c_bool)
    c%converge_on_zero_diff = logical(f%converge_on_zero_diff, c_bool)
    c%fcn_count = f%fcn_count
    c%gradient_count = f%gradient_count
    c%iter_count = f%iter_count
    c%jacobian_count = f%jacobian_count
end subroutine


pure subroutine convert_from_c_iteration_behavior(f, c)
    type(iteration_behavior), intent(out) :: f
    type(c_iteration_behavior), intent(in) :: c
    f%converge_on_chng = c%converge_on_chng
    f%converge_on_fcn = c%converge_on_fcn
    f%converge_on_zero_diff = c%converge_on_zero_diff
    f%fcn_count = c%fcn_count
    f%gradient_count = c%gradient_count
    f%iter_count = c%iter_count
    f%jacobian_count = c%jacobian_count
end subroutine


pure function cw_eval(this, bin) result(rst)
    class(c_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    rst = this%fcn(this%size, bin)
end function


pure subroutine convert_to_c_quaternion(qc, qf)
    type(c_quaternion), intent(out) :: qc
    type(quaternion), intent(in) :: qf
    qc%w = qf%w
    qc%x = qf%x
    qc%y = qf%y
    qc%z = qf%z
end subroutine


pure subroutine convert_from_c_quaternion(qf, qc)
    type(quaternion), intent(out) :: qf
    type(c_quaternion), intent(in) :: qc
    qf%w = qc%w
    qf%x = qc%x
    qf%y = qc%y
    qf%z = qc%z
end subroutine


pure subroutine convert_to_c_line(lc, lf)
    type(c_line), intent(out) :: lc
    type(line), intent(in) :: lf
    lc%r0 = lf%r0
    lc%v = lf%v
end subroutine


pure subroutine convert_from_c_line(lf, lc)
    type(line), intent(out) :: lf
    type(c_line), intent(in) :: lc
    lf%r0 = lc%r0
    lf%v = lc%v
end subroutine


pure subroutine convert_to_c_plane(pc, pf)
    type(c_plane), intent(out) :: pc
    type(plane), intent(in) :: pf
    pc%a = pf%a
    pc%b = pf%b
    pc%c = pf%c
    pc%d = pf%d
end subroutine


pure subroutine convert_from_c_plane(pf, pc)
    type(plane), intent(out) :: pf
    type(c_plane), intent(in) :: pc
    pf%a = pc%a
    pf%b = pc%b
    pf%c = pc%c
    pf%d = pc%d
end subroutine


pure subroutine convert_to_c_plucker_line(pc, pf)
    type(c_plucker_line), intent(out) :: pc
    type(plucker_line), intent(in) :: pf
    pc%u = pf%v(1:3)
    pc%m = pf%v(4:6)
end subroutine


pure subroutine convert_from_c_plucker_line(pf, pc)
    type(plucker_line), intent(out) :: pf
    type(c_plucker_line), intent(in) :: pc
    pf%v(1:3) = pc%u
    pf%v(4:6) = pc%m
end subroutine


pure subroutine convert_to_c_coordinate_system(cc, cf)
    type(c_coordinate_system), intent(out) :: cc
    type(coordinate_system), intent(in) :: cf
    cc%origin = cf%origin
    cc%i = cf%i
    cc%j = cf%j
    cc%k = cf%k
end subroutine


pure subroutine convert_from_c_coordinate_system(cf, cc)
    type(coordinate_system), intent(out) :: cf
    type(c_coordinate_system), intent(in) :: cc
    cf%origin = cc%origin
    cf%i = cc%i
    cf%j = cc%j
    cf%k = cc%k
end subroutine


pure subroutine convert_to_c_dh_parameter_set(dc, df)
    type(c_dh_parameter_set), intent(out) :: dc
    type(dh_parameter_set), intent(in) :: df
    dc%joint_angle = df%joint_angle
    dc%link_length = df%link_length
    dc%link_offset = df%link_offset
    dc%link_twist = df%link_twist
end subroutine


pure subroutine convert_from_c_dh_parameter_set(df, dc)
    type(dh_parameter_set), intent(out) :: df
    type(c_dh_parameter_set), intent(in) :: dc
    df%joint_angle = dc%joint_angle
    df%link_length = dc%link_length
    df%link_offset = dc%link_offset
    df%link_twist = dc%link_twist
end subroutine


subroutine convert_from_c_dh_table(df, dc)
    type(dh_table), intent(out) :: df
    type(c_dh_table), intent(in) :: dc
    integer(int32) :: i, n
    type(c_dh_parameter_set), pointer, dimension(:) :: ptr
    n = int(dc%count, int32)
    allocate(df%parameters(n))
    call c_f_pointer(dc%parameters, ptr, [n])
    do i = 1, n
        df%parameters(i) = ptr(i)   ! = automatically makes the conversion
    end do
    ! Notes: No memory need be allocated for ptr as it is just a container for
    ! the already allocated C pointer.
end subroutine


subroutine convert_to_c_dh_table(dc, df) ! cannot be used by assignment operator
    type(c_dh_table), intent(inout) :: dc
    type(dh_table), intent(in) :: df
    integer(int32) :: i, n
    type(c_dh_parameter_set), pointer, dimension(:) :: ptr
    n = size(df%parameters)
    dc%count = int(n, c_int)
    call c_f_pointer(dc%parameters, ptr, [n])
    do i = 1, n
        ptr(i) = df%parameters(i)
    end do
end subroutine


pure subroutine convert_to_c_binary_link(bc, bf)
    type(c_binary_link), intent(out) :: bc
    type(binary_link), intent(in) :: bf
    bc%joint_angle = bf%joint_angle
    bc%joint_type = bf%joint_type
    bc%link_length = bf%link_length
    bc%link_offset = bf%link_offset
    bc%link_twist = bf%link_twist
    bc%mass = bf%mass
    bc%cg = bf%cg
    bc%inertia = reshape(bf%inertia, [9])
end subroutine


pure subroutine convert_from_c_binary_link(bf, bc)
    type(binary_link), intent(out) :: bf
    type(c_binary_link), intent(in) :: bc
    bf%joint_angle = bc%joint_angle
    bf%joint_type = bc%joint_type
    bf%link_length = bc%link_length
    bf%link_offset = bc%link_offset
    bf%link_twist = bc%link_twist
    bf%mass = bc%mass
    bf%cg = bc%cg
    bf%inertia = reshape(bc%inertia, [3, 3])
end subroutine


subroutine convert_from_c_serial_linkage(sf, sc)
    type(serial_linkage), intent(out) :: sf
    type(c_serial_linkage), intent(in) :: sc
    integer(int32) :: i, n
    type(binary_link), allocatable, dimension(:) :: links
    type(c_binary_link), pointer, dimension(:) :: ptr
    n = int(sc%link_count, int32)
    allocate(links(n))
    call c_f_pointer(sc%links, ptr, [n])
    do i = 1, n
        links(i) = ptr(i)
    end do
    sf = serial_linkage(links)
end subroutine


subroutine convert_to_c_serial_linkage(sc, sf)
    type(c_serial_linkage), intent(out) :: sc
    type(serial_linkage), intent(in) :: sf
    integer(int32) :: i, n
    type(c_binary_link), pointer, dimension(:) :: ptr
    n = sf%get_link_count()
    sc%link_count = int(n, c_int)
    call c_f_pointer(sc%links, ptr, [n])
    do i = 1, n
        ptr(i) = sf%get_link(i)
    end do
end subroutine


subroutine convert_to_c_polynomial(pc, pf)
    type(c_polynomial), intent(out) :: pc
    type(polynomial), intent(in) :: pf

    integer(c_int) :: i, flag, n
    real(c_double), pointer, dimension(:) :: ptr

    flag = c_alloc_polynomial(pf%order(), pc)
    if (flag < 0) return
    n = pf%order() + 1
    call c_f_pointer(pc%coefficients, ptr, [n])
    do i = 1, n
        ptr(i) = pf%get(i)
    end do
end subroutine


subroutine convert_from_c_polynomial(pf, pc)
    type(polynomial), intent(out) :: pf
    type(c_polynomial), intent(in) :: pc

    integer(c_int) :: i, order, n
    real(c_double), pointer, dimension(:) :: ptr

    order = pc%order
    n = order + 1
    call c_f_pointer(pc%coefficients, ptr, [n])
    call pf%initialize(ptr)
end subroutine


subroutine convert_to_c_transfer_function(tc, tf)
    type(c_transfer_function), intent(out) :: tc
    type(transfer_function), intent(in) :: tf

    tc%numerator = tf%Y
    tc%denominator = tf%X
end subroutine


subroutine convert_from_c_transfer_function(tf, tc)
    type(transfer_function), intent(out) :: tf
    type(c_transfer_function), intent(in) :: tc

    tf%Y = tc%numerator
    tf%X = tc%denominator
end subroutine


subroutine convert_to_c_state_space(sc, sf)
    type(c_state_space_model), intent(out) :: sc
    type(state_space), intent(in) :: sf

    integer(int32) :: i, j, n, m, p, flag
    real(real64), pointer, dimension(:,:) :: a, b, c, d

    n = size(sf%A, 1)   ! size
    m = size(sf%B, 2)   ! # of inputs
    p = size(sf%C, 1)   ! # of outputs
    flag = c_alloc_state_space_model(n, m, p, sc)
    if (flag < 0) return
    call c_f_pointer(sc%A, a, [n, n])
    call c_f_pointer(sc%B, b, [n, m])
    call c_f_pointer(sc%C, c, [p, n])
    call c_f_pointer(sc%D, d, [p, m])
    do j = 1, n
        do i = 1, n
            a(i,j) = sf%A(i,j)
        end do
        do i = 1, p
            c(i,j) = sf%C(i,j)
        end do
    end do
    do j = 1, m
        do i = 1, n
            b(i,j) = sf%B(i,j)
        end do
        do i = 1, p
            d(i,j) = sf%D(i,j)
        end do
    end do
end subroutine


subroutine convert_from_c_state_space(sf, sc)
    type(state_space), intent(out) :: sf
    type(c_state_space_model), intent(in) :: sc

    integer(int32) :: n, m, p
    real(real64), pointer, dimension(:,:) :: a, b, c, d

    n = sc%dimension
    m = sc%n_inputs
    p = sc%n_outputs
    call c_f_pointer(sc%A, a, [n, n])
    call c_f_pointer(sc%B, b, [n, m])
    call c_f_pointer(sc%C, c, [p, n])
    call c_f_pointer(sc%D, d, [p, m])
    sf = state_space(a, b, c, d)
end subroutine


end module
