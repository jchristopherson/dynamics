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
module dynamics_flexible_discrete
    !! Two-node spring and viscous damper elements in two and three dimensions.
    !! The spring uses a nonlinear length-based energy.  The damper uses a
    !! quadratic dissipation function of the relative nodal velocity along the
    !! current element axis.
    use iso_fortran_env
    use autodiff
    use dynamics_flexible_types
    implicit none
    private
    public :: flexible_spring_2d
    public :: flexible_spring_3d
    public :: flexible_damper_2d
    public :: flexible_damper_3d

    type, extends(finite_element) :: flexible_spring_2d
        !! A two-dimensional spring with energy
        !! \(W(\delta) = \frac{1}{2} k \delta^2 + \frac{1}{4} k_3 \delta^4\),
        !! where \(\delta = L - L_0\) is the change in length.
        real(real64) :: reference_coordinates(2,2) = 0.0d0
            !! The reference nodal coordinates, one column per node.
        real(real64) :: stiffness = 1.0d0
            !! The linear stiffness, \(k\).
        real(real64) :: cubic_stiffness = 0.0d0
            !! The nonlinear stiffness, \(k_3\).
    contains
        procedure :: get_dof_count => spring_2d_dof_count
        procedure :: strain_energy => spring_2d_strain_energy
    end type

    type, extends(finite_element) :: flexible_spring_3d
        !! A three-dimensional spring with energy
        !! \(W(\delta) = \frac{1}{2} k \delta^2 + \frac{1}{4} k_3 \delta^4\),
        !! where \(\delta = L - L_0\) is the change in length.
        real(real64) :: reference_coordinates(3,2) = 0.0d0
            !! The reference nodal coordinates, one column per node.
        real(real64) :: stiffness = 1.0d0
            !! The linear stiffness, \(k\).
        real(real64) :: cubic_stiffness = 0.0d0
            !! The nonlinear stiffness, \(k_3\).
    contains
        procedure :: get_dof_count => spring_3d_dof_count
        procedure :: strain_energy => spring_3d_strain_energy
    end type

    type, extends(dissipative_element) :: flexible_damper_2d
        !! A two-dimensional viscous damper with dissipation
        !! \(D = \frac{1}{2} c \dot{L}^2\).
        real(real64) :: reference_coordinates(2,2) = 0.0d0
            !! The reference nodal coordinates, one column per node.
        real(real64) :: damping = 1.0d0
            !! The viscous damping coefficient, \(c\).
    contains
        procedure :: get_dof_count => damper_2d_dof_count
        procedure :: dissipation => damper_2d_dissipation
    end type

    type, extends(dissipative_element) :: flexible_damper_3d
        !! A three-dimensional viscous damper with dissipation
        !! \(D = \frac{1}{2} c \dot{L}^2\).
        real(real64) :: reference_coordinates(3,2) = 0.0d0
            !! The reference nodal coordinates, one column per node.
        real(real64) :: damping = 1.0d0
            !! The viscous damping coefficient, \(c\).
    contains
        procedure :: get_dof_count => damper_3d_dof_count
        procedure :: dissipation => damper_3d_dissipation
    end type

contains
! ------------------------------------------------------------------------------
function spring_2d_dof_count(this) result(n)
    !! Returns the number of degrees of freedom of the 2D spring, which is four.
    class(flexible_spring_2d), intent(in) :: this
        !! The element.
    integer(int32) :: n
        !! The number of element degrees of freedom.

    n = 4
end function

! ------------------------------------------------------------------------------
function spring_3d_dof_count(this) result(n)
    !! Returns the number of degrees of freedom of the 3D spring, which is six.
    class(flexible_spring_3d), intent(in) :: this
        !! The element.
    integer(int32) :: n
        !! The number of element degrees of freedom.

    n = 6
end function

! ------------------------------------------------------------------------------
function damper_2d_dof_count(this) result(n)
    !! Returns the number of degrees of freedom of the 2D damper, which is four.
    class(flexible_damper_2d), intent(in) :: this
        !! The element.
    integer(int32) :: n
        !! The number of element degrees of freedom.

    n = 4
end function

! ------------------------------------------------------------------------------
function damper_3d_dof_count(this) result(n)
    !! Returns the number of degrees of freedom of the 3D damper, which is six.
    class(flexible_damper_3d), intent(in) :: this
        !! The element.
    integer(int32) :: n
        !! The number of element degrees of freedom.

    n = 6
end function

! ------------------------------------------------------------------------------
function spring_2d_strain_energy(this, u) result(w)
    !! Evaluates the strain energy of the 2D spring for hyper-dual displacements.
    class(flexible_spring_2d), intent(in) :: this
        !! The element.
    type(hyperdual), intent(in), dimension(:) :: u
        !! The hyper-dual nodal displacements, ordered by node.
    type(hyperdual) :: w
        !! The strain energy.

    w = spring_energy(this%reference_coordinates, this%stiffness, &
        this%cubic_stiffness, u)
end function

! ------------------------------------------------------------------------------
function spring_3d_strain_energy(this, u) result(w)
    !! Evaluates the strain energy of the 3D spring for hyper-dual displacements.
    class(flexible_spring_3d), intent(in) :: this
        !! The element.
    type(hyperdual), intent(in), dimension(:) :: u
        !! The hyper-dual nodal displacements, ordered by node.
    type(hyperdual) :: w
        !! The strain energy.

    w = spring_energy(this%reference_coordinates, this%stiffness, &
        this%cubic_stiffness, u)
end function

! ------------------------------------------------------------------------------
function damper_2d_dissipation(this, u, v) result(d)
    !! Evaluates the dissipation of the 2D damper for hyper-dual velocities.
    class(flexible_damper_2d), intent(in) :: this
        !! The element.
    real(real64), intent(in), dimension(:) :: u
        !! The nodal displacements, which set the current axis direction.
    type(hyperdual), intent(in), dimension(:) :: v
        !! The hyper-dual nodal velocities, ordered by node.
    type(hyperdual) :: d
        !! The dissipation.

    d = damper_dissipation(this%reference_coordinates, this%damping, u, v)
end function

! ------------------------------------------------------------------------------
function damper_3d_dissipation(this, u, v) result(d)
    !! Evaluates the dissipation of the 3D damper for hyper-dual velocities.
    class(flexible_damper_3d), intent(in) :: this
        !! The element.
    real(real64), intent(in), dimension(:) :: u
        !! The nodal displacements, which set the current axis direction.
    type(hyperdual), intent(in), dimension(:) :: v
        !! The hyper-dual nodal velocities, ordered by node.
    type(hyperdual) :: d
        !! The dissipation.

    d = damper_dissipation(this%reference_coordinates, this%damping, u, v)
end function

! ------------------------------------------------------------------------------
function spring_energy(x0, k, k3, u) result(w)
    !! Strain energy of a spring with reference coordinates `x0`.  The current
    !! length is computed from the hyper-dual displacements, and the change in
    !! length \(\delta = L - L_0\) gives \(W = \frac{1}{2} k \delta^2 +
    !! \frac{1}{4} k_3 \delta^4\).
    real(real64), intent(in), dimension(:,:) :: x0
        !! The reference nodal coordinates, one column per node.
    real(real64), intent(in) :: k
        !! The linear stiffness.
    real(real64), intent(in) :: k3
        !! The nonlinear stiffness.
    type(hyperdual), intent(in), dimension(:) :: u
        !! The hyper-dual nodal displacements.
    type(hyperdual) :: w
        !! The strain energy.

    integer(int32) :: i, d
    type(hyperdual) :: l2, len, delta, dx

    d = size(x0, 1)
    l2 = 0.0d0
    do i = 1, d
        dx = (x0(i,2) - x0(i,1)) + (u(d + i) - u(i))
        l2 = l2 + dx * dx
    end do
    len = sqrt(l2)
    delta = len - sqrt(sum((x0(:,2) - x0(:,1))**2))
    w = 0.5d0 * k * delta * delta + 0.25d0 * k3 * delta**4
end function

! ------------------------------------------------------------------------------
function damper_dissipation(x0, c, u, v) result(d)
    !! Dissipation of a viscous damper, \(D = \frac{1}{2} c \dot{L}^2\).  The
    !! current axis direction is taken from the real displacements `u`, and the
    !! rate of change of length is the projection of the hyper-dual relative
    !! velocity onto that direction.
    real(real64), intent(in), dimension(:,:) :: x0
        !! The reference nodal coordinates, one column per node.
    real(real64), intent(in) :: c
        !! The viscous damping coefficient.
    real(real64), intent(in), dimension(:) :: u
        !! The nodal displacements, which set the current axis direction.
    type(hyperdual), intent(in), dimension(:) :: v
        !! The hyper-dual nodal velocities.
    type(hyperdual) :: d
        !! The dissipation.

    integer(int32) :: i, n
    real(real64) :: dir(size(x0, 1)), length
    type(hyperdual) :: rate

    n = size(x0, 1)
    dir = x0(:,2) - x0(:,1) + u(n+1:2*n) - u(1:n)
    length = sqrt(sum(dir**2))
    rate = 0.0d0
    do i = 1, n
        rate = rate + dir(i) * (v(n + i) - v(i))
    end do
    rate = rate / length
    d = 0.5d0 * c * rate * rate
end function

! ------------------------------------------------------------------------------
end module
