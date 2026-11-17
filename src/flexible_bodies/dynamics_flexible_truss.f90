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
module dynamics_flexible_truss
    !! Two-node truss elements in two and three dimensions with a finite-strain
    !! (Green-Lagrange) axial strain measure.
    use iso_fortran_env
    use autodiff
    use dynamics_flexible_types
    use dynamics_flexible_materials
    implicit none
    private
    public :: flexible_truss_2d
    public :: flexible_truss_3d

    type, extends(finite_element) :: flexible_truss_2d
        !! A two-dimensional truss element with two nodes and four degrees of
        !! freedom ordered as \([u_{1x}, u_{1y}, u_{2x}, u_{2y}]\).
        real(real64) :: reference_coordinates(2,2) = 0.0d0
            !! The reference nodal coordinates, one column per node.
        real(real64) :: area = 1.0d0
            !! The cross-sectional area.
        class(flexible_material), allocatable :: material
            !! The material model.
    contains
        procedure :: get_dof_count => truss_2d_dof_count
        procedure :: strain_energy => truss_2d_strain_energy
    end type

    type, extends(finite_element) :: flexible_truss_3d
        !! A three-dimensional truss element with two nodes and six degrees of
        !! freedom ordered as \([u_{1x}, u_{1y}, u_{1z}, u_{2x}, u_{2y}, u_{2z}]\).
        real(real64) :: reference_coordinates(3,2) = 0.0d0
            !! The reference nodal coordinates, one column per node.
        real(real64) :: area = 1.0d0
            !! The cross-sectional area.
        class(flexible_material), allocatable :: material
            !! The material model.
    contains
        procedure :: get_dof_count => truss_3d_dof_count
        procedure :: strain_energy => truss_3d_strain_energy
    end type

contains
! ------------------------------------------------------------------------------
function truss_2d_dof_count(this) result(n)
    class(flexible_truss_2d), intent(in) :: this
    integer(int32) :: n

    n = 4
end function

! ------------------------------------------------------------------------------
function truss_3d_dof_count(this) result(n)
    class(flexible_truss_3d), intent(in) :: this
    integer(int32) :: n

    n = 6
end function

! ------------------------------------------------------------------------------
function truss_2d_strain_energy(this, u) result(w)
    class(flexible_truss_2d), intent(in) :: this
    type(hyperdual), intent(in), dimension(:) :: u
    type(hyperdual) :: w

    w = bar_strain_energy(this%reference_coordinates, this%area, this%material, u)
end function

! ------------------------------------------------------------------------------
function truss_3d_strain_energy(this, u) result(w)
    class(flexible_truss_3d), intent(in) :: this
    type(hyperdual), intent(in), dimension(:) :: u
    type(hyperdual) :: w

    w = bar_strain_energy(this%reference_coordinates, this%area, this%material, u)
end function

! ------------------------------------------------------------------------------
function bar_strain_energy(x0, area, mat, u) result(w)
    !! Strain energy of a truss with reference coordinates `x0` (one column per
    !! node).  The Green-Lagrange strain is \(E = (L^2 - L_0^2) / (2 L_0^2)\).
    real(real64), intent(in), dimension(:,:) :: x0
    real(real64), intent(in) :: area
    class(flexible_material), intent(in) :: mat
    type(hyperdual), intent(in), dimension(:) :: u
    type(hyperdual) :: w

    integer(int32) :: i, d
    real(real64) :: l0sq
    type(hyperdual) :: l2, e, dx

    d = size(x0, 1)
    l0sq = sum((x0(:,2) - x0(:,1))**2)
    l2 = 0.0d0
    do i = 1, d
        dx = (x0(i,2) - x0(i,1)) + (u(d + i) - u(i))
        l2 = l2 + dx * dx
    end do
    e = (l2 - l0sq) / (2.0d0 * l0sq)
    w = sqrt(l0sq) * area * mat%strain_energy_density(e)
end function

! ------------------------------------------------------------------------------
end module
