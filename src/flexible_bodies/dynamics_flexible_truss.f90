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
    !! The element integrates the material's axial energy density over the
    !! reference volume; it has no bending or shear degrees of freedom.
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_is_finite
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
            !! Reference nodal coordinates, one column per node; the nodes must be distinct.
        real(real64) :: area = 1.0d0
            !! Reference cross-sectional area, in units of length squared.
        class(flexible_material), allocatable :: material
            !! Axial material law; must be allocated before evaluating the energy.
    contains
        procedure :: get_dof_count => truss_2d_dof_count
        procedure :: strain_energy => truss_2d_strain_energy
        procedure :: mass_matrix => truss_2d_mass_matrix
    end type

    type, extends(finite_element) :: flexible_truss_3d
        !! A three-dimensional truss element with two nodes and six degrees of
        !! freedom ordered as \([u_{1x}, u_{1y}, u_{1z}, u_{2x}, u_{2y}, u_{2z}]\).
        real(real64) :: reference_coordinates(3,2) = 0.0d0
            !! Reference nodal coordinates, one column per node; the nodes must be distinct.
        real(real64) :: area = 1.0d0
            !! Reference cross-sectional area, in units of length squared.
        class(flexible_material), allocatable :: material
            !! Axial material law; must be allocated before evaluating the energy.
    contains
        procedure :: get_dof_count => truss_3d_dof_count
        procedure :: strain_energy => truss_3d_strain_energy
        procedure :: mass_matrix => truss_3d_mass_matrix
    end type

contains
function truss_2d_mass_matrix(this) result(m)
    !! Consistent translational mass for a two-node planar truss.
    class(flexible_truss_2d), intent(in) :: this
    real(real64), allocatable :: m(:,:)

    m = bar_mass_matrix(this%reference_coordinates, this%area, this%material)
end function

function truss_3d_mass_matrix(this) result(m)
    !! Consistent translational mass for a two-node spatial truss.
    class(flexible_truss_3d), intent(in) :: this
    real(real64), allocatable :: m(:,:)

    m = bar_mass_matrix(this%reference_coordinates, this%area, this%material)
end function

function bar_mass_matrix(x, area, material) result(m)
    !! Integrates the product of linear nodal shape functions over the reference bar.
    real(real64), intent(in) :: x(:,:), area
    class(flexible_material), allocatable, intent(in) :: material
    real(real64) :: m(2*size(x,1),2*size(x,1)), length, mass_sixth
    integer :: i, n

    if (.not. allocated(material)) error stop "Flexible truss requires a material"
    length = norm2(x(:,2) - x(:,1))
    if (.not. ieee_is_finite(length) .or. length <= 0.0d0 .or. &
        .not. ieee_is_finite(area) .or. area <= 0.0d0) error stop "Invalid flexible truss geometry"
    if (.not. ieee_is_finite(material%density) .or. material%density < 0.0d0) &
        error stop "Invalid flexible truss density"
    mass_sixth = material%density * area * length / 6.0d0
    n = size(x,1)
    m = 0.0d0
    do i = 1, n
        m(i,i) = 2.0d0 * mass_sixth
        m(n+i,n+i) = m(i,i)
        m(i,n+i) = mass_sixth
        m(n+i,i) = mass_sixth
    end do
end function

! ------------------------------------------------------------------------------
function truss_2d_dof_count(this) result(n)
    !! Returns the number of degrees of freedom of the 2D truss, which is four.
    class(flexible_truss_2d), intent(in) :: this
        !! The element.
    integer(int32) :: n
        !! The number of element degrees of freedom.

    n = 4
end function

! ------------------------------------------------------------------------------
function truss_3d_dof_count(this) result(n)
    !! Returns the number of degrees of freedom of the 3D truss, which is six.
    class(flexible_truss_3d), intent(in) :: this
        !! The element.
    integer(int32) :: n
        !! The number of element degrees of freedom.

    n = 6
end function

! ------------------------------------------------------------------------------
function truss_2d_strain_energy(this, u) result(w)
    !! Evaluates the strain energy of the 2D truss for hyper-dual displacements.
    class(flexible_truss_2d), intent(in) :: this
        !! The element.
    type(hyperdual), intent(in), dimension(:) :: u
        !! Four hyper-dual translational DOFs, ordered by node.
    type(hyperdual) :: w
        !! The strain energy.

    w = bar_strain_energy(this%reference_coordinates, this%area, this%material, u)
end function

! ------------------------------------------------------------------------------
function truss_3d_strain_energy(this, u) result(w)
    !! Evaluates the strain energy of the 3D truss for hyper-dual displacements.
    class(flexible_truss_3d), intent(in) :: this
        !! The element.
    type(hyperdual), intent(in), dimension(:) :: u
        !! Six hyper-dual translational DOFs, ordered by node.
    type(hyperdual) :: w
        !! The strain energy.

    w = bar_strain_energy(this%reference_coordinates, this%area, this%material, u)
end function

! ------------------------------------------------------------------------------
function bar_strain_energy(x0, area, mat, u) result(w)
    !! Strain energy of a truss with reference coordinates `x0` (one column per
    !! node).  The Green-Lagrange strain is \(E = (L^2 - L_0^2) / (2 L_0^2)\),
    !! and the energy is \(W = L_0 A W_m(E)\), where \(W_m\) is the material
    !! strain energy per unit reference volume. The reference length \(L_0\)
    !! must be nonzero.
    real(real64), intent(in), dimension(:,:) :: x0
        !! The reference nodal coordinates, one column per node.
    real(real64), intent(in) :: area
        !! Reference cross-sectional area.
    class(flexible_material), intent(in) :: mat
        !! Material supplying the axial law for the Green-Lagrange strain.
    type(hyperdual), intent(in), dimension(:) :: u
        !! Hyper-dual translational DOFs, with `size(x0,1)` components per node.
    type(hyperdual) :: w
        !! The strain energy.

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
