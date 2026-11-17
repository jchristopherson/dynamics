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
    contains
        procedure(dof_count_function), deferred :: get_dof_count
            !! Returns the number of element degrees of freedom.
    end type

    type, extends(flexible_element), abstract :: finite_element
        !! Base type for elements that store strain energy.
    contains
        procedure(strain_energy_function), deferred :: strain_energy
            !! Evaluates the strain energy for hyper-dual element displacements.
        procedure :: internal_force => finite_element_internal_force
            !! Computes the internal force vector and tangent stiffness matrix.
    end type

    type, extends(flexible_element), abstract :: dissipative_element
        !! Base type for elements that store a dissipation function.
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
        real(real64), allocatable, dimension(:) :: u
    end type

    abstract interface
        function dof_count_function(this) result(n)
            import
            class(flexible_element), intent(in) :: this
            integer(int32) :: n
        end function

        function strain_energy_function(this, u) result(w)
            import
            class(finite_element), intent(in) :: this
            type(hyperdual), intent(in), dimension(:) :: u
                !! Hyper-dual element displacements.
            type(hyperdual) :: w
        end function

        function dissipation_function(this, u, v) result(d)
            import
            class(dissipative_element), intent(in) :: this
            real(real64), intent(in), dimension(:) :: u
                !! Element displacements.
            type(hyperdual), intent(in), dimension(:) :: v
                !! Hyper-dual element velocities.
            type(hyperdual) :: d
        end function
    end interface

contains
! ------------------------------------------------------------------------------
subroutine finite_element_internal_force(this, u, f, k)
    !! Computes the internal force vector `f` and tangent stiffness matrix `k`
    !! at the element displacements `u`.
    class(finite_element), intent(in) :: this
    real(real64), intent(in), dimension(:) :: u
    real(real64), intent(out), dimension(:) :: f
    real(real64), intent(out), dimension(:,:) :: k

    call hessian(strain_energy_wrapper, u, k, g = f, args = this)
end subroutine

! ------------------------------------------------------------------------------
function strain_energy_wrapper(x, args) result(w)
    !! Adapts the strain energy method to the hyper-dual function interface.
    type(hyperdual), intent(in), dimension(:) :: x
    class(*), intent(in), optional :: args
    type(hyperdual) :: w

    select type (args)
    class is (finite_element)
        w = args%strain_energy(x)
    class default
        w = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
subroutine dissipative_element_damping_matrix(this, u, c)
    !! Computes the damping matrix `c` at the element displacements `u`.  The
    !! dissipation is quadratic in the velocities, so the Hessian is independent
    !! of the evaluation point in velocity space.
    class(dissipative_element), intent(in) :: this
    real(real64), intent(in), dimension(:) :: u
    real(real64), intent(out), dimension(:,:) :: c

    type(dissipation_context) :: ctx
    real(real64), allocatable, dimension(:) :: v

    allocate(ctx%element, source = this)
    ctx%u = u
    allocate(v(size(u)), source = 0.0d0)
    call hessian(dissipation_wrapper, v, c, args = ctx)
end subroutine

! ------------------------------------------------------------------------------
function dissipation_wrapper(x, args) result(d)
    !! Adapts the dissipation method to the hyper-dual function interface.
    type(hyperdual), intent(in), dimension(:) :: x
    class(*), intent(in), optional :: args
    type(hyperdual) :: d

    select type (args)
    type is (dissipation_context)
        d = args%element%dissipation(args%u, x)
    class default
        d = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
end module
