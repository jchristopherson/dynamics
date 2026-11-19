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
module dynamics_flexible_types
    !! Base types for flexible-body elements.  A finite element supplies its
    !! strain energy as a function of hyper-dual nodal displacements.  The
    !! internal force vector and tangent stiffness matrix are obtained from the
    !! energy by automatic differentiation.  A dissipative element supplies a
    !! dissipation function of the nodal velocities in the same manner, and its
    !! damping matrix is obtained by automatic differentiation.
    use iso_fortran_env
    use autodiff
    implicit none
    private
    public :: flexible_element
    public :: finite_element
    public :: dissipative_element
    public :: dissipation_context

    type, abstract :: flexible_element
        !! Base type for all flexible-body elements.
        !! Concrete elements define the ordering of their degrees of freedom;
        !! all input and output arrays use that same ordering.
    contains
        procedure(dof_count_function), deferred :: get_dof_count
            !! Returns the number of element degrees of freedom.
        procedure :: mass_matrix => massless_element_matrix
            !! Returns the reference mass matrix in element DOF order. The
            !! default is zero for massless connectors and legacy elements.
            !! An override must return an allocated, symmetric square matrix
            !! with one row and column per element DOF.
    end type

    type, extends(flexible_element), abstract :: finite_element
        !! Base type for elements that store strain energy.
        !! The internal force uses the positive energy-gradient convention.
    contains
        procedure(strain_energy_function), deferred :: strain_energy
            !! Evaluates the strain energy for hyper-dual element displacements.
        procedure :: internal_force => finite_element_internal_force
            !! Computes the internal force vector and tangent stiffness matrix.
    end type

    type, extends(flexible_element), abstract :: dissipative_element
        !! Base type for elements that store a dissipation function.
        !! The default damping matrix is the velocity Hessian at zero velocity.
    contains
        procedure(dissipation_function), deferred :: dissipation
            !! Evaluates the dissipation for displacements and hyper-dual velocities.
        procedure :: damping_matrix => dissipative_element_damping_matrix
            !! Computes the damping matrix from the dissipation function.
    end type

    type :: dissipation_context
        !! Caller-owned context passed to the dissipation function during
        !! automatic differentiation.
        class(dissipative_element), allocatable :: element
            !! Allocated copy of the element supplying the dissipation.
        real(real64), allocatable, dimension(:) :: u
            !! Displacements held fixed while differentiating in velocity.
    end type

    abstract interface
        function dof_count_function(this) result(n)
            !! Returns the number of degrees of freedom of the element.
            import
            class(flexible_element), intent(in) :: this
                !! The element.
            integer(int32) :: n
                !! The number of element degrees of freedom.
        end function

        function strain_energy_function(this, u) result(w)
            !! Evaluates the strain energy of the element.  The function must
            !! be written for hyper-dual arguments so that its derivatives can
            !! be obtained by automatic differentiation.
            import
            class(finite_element), intent(in) :: this
                !! The element.
            type(hyperdual), intent(in), dimension(:) :: u
                !! Hyper-dual generalized displacements in element DOF order.
            type(hyperdual) :: w
                !! Total stored strain energy of the element.
        end function

        function dissipation_function(this, u, v) result(d)
            !! Evaluates the dissipation function of the element.  The function
            !! must be written for hyper-dual velocities so that the damping
            !! matrix can be obtained by automatic differentiation.
            import
            class(dissipative_element), intent(in) :: this
                !! The element.
            real(real64), intent(in), dimension(:) :: u
                !! Element displacements, held at their evaluation point.
            type(hyperdual), intent(in), dimension(:) :: v
                !! Hyper-dual generalized velocities in the same DOF order as `u`.
            type(hyperdual) :: d
                !! Element dissipation potential; its velocity gradient is the damping force.
        end function
    end interface

contains
function massless_element_matrix(this) result(m)
    !! Returns a zero mass matrix; override for elements that carry inertia.
    class(flexible_element), intent(in) :: this
        !! Element supplying the size and ordering of the matrix.
    real(real64), allocatable, dimension(:,:) :: m
        !! Zero matrix with one row and column per element DOF.
    integer :: n

    n = this%get_dof_count()
    allocate(m(n,n), source=0.0d0)
end function

! ------------------------------------------------------------------------------
subroutine finite_element_internal_force(this, u, f, k)
    !! Computes the internal force vector and tangent stiffness matrix of the
    !! element at the displacements `u`.  The force is the gradient of the
    !! strain energy and the tangent is its Hessian, both obtained by automatic
    !! differentiation.  The tangent costs \(n(n+1)/2\) energy evaluations for
    !! an element with \(n\) degrees of freedom.
    class(finite_element), intent(in) :: this
        !! The element.
    real(real64), intent(in), dimension(:) :: u
        !! The element displacements, of length `get_dof_count()`.
    real(real64), intent(out), dimension(:) :: f
        !! Energy gradient, of length `get_dof_count()`.
    real(real64), intent(out), dimension(:,:) :: k
        !! Energy Hessian, an \(n \times n\) matrix in element DOF order.

    call hessian(strain_energy_wrapper, u, k, g = f, args = this)
end subroutine

! ------------------------------------------------------------------------------
function strain_energy_wrapper(x, args) result(w)
    !! Adapts the `strain_energy` type-bound method to the hyper-dual function
    !! interface expected by the `hessian` driver.  The element is passed
    !! through `args`.  If `args` does not hold a `finite_element`, the result
    !! is zero.
    type(hyperdual), intent(in), dimension(:) :: x
        !! The hyper-dual element displacements.
    class(*), intent(in), optional :: args
        !! The finite element that supplies the strain energy.
    type(hyperdual) :: w
        !! The strain energy.

    select type (args)
    class is (finite_element)
        w = args%strain_energy(x)
    class default
        w = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
subroutine dissipative_element_damping_matrix(this, u, c)
    !! Computes the damping matrix of the element at the displacements `u`.
    !! The Hessian of the dissipation with respect to the velocities is taken
    !! by automatic differentiation.  The dissipation is assumed quadratic in
    !! the velocities, so the result does not depend on the velocity at which
    !! it is evaluated, and the velocity is set to zero.
    class(dissipative_element), intent(in) :: this
        !! The element.
    real(real64), intent(in), dimension(:) :: u
        !! The element displacements, of length `get_dof_count()`.
    real(real64), intent(out), dimension(:,:) :: c
        !! Velocity Hessian, an \(n \times n\) matrix with \(n\) equal to `get_dof_count()`.

    type(dissipation_context) :: ctx
    real(real64), allocatable, dimension(:) :: v

    allocate(ctx%element, source = this)
    ctx%u = u
    allocate(v(size(u)), source = 0.0d0)
    call hessian(dissipation_wrapper, v, c, args = ctx)
end subroutine

! ------------------------------------------------------------------------------
function dissipation_wrapper(x, args) result(d)
    !! Adapts the `dissipation` type-bound method to the hyper-dual function
    !! interface expected by the `hessian` driver.  The element and the real
    !! displacements are taken from the `dissipation_context` in `args`.  If
    !! `args` does not hold a `dissipation_context`, the result is zero.
    type(hyperdual), intent(in), dimension(:) :: x
        !! The hyper-dual element velocities.
    class(*), intent(in), optional :: args
        !! Dissipation context with an allocated element and fixed displacements.
    type(hyperdual) :: d
        !! The dissipation.

    select type (args)
    type is (dissipation_context)
        d = args%element%dissipation(args%u, x)
    class default
        d = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
end module
