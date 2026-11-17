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
module dynamics_flexible_materials
    !! Material models for flexible-body elements.  A material supplies its
    !! strain energy density as a function of a Green-Lagrange strain, written
    !! for hyper-dual inputs so that the element stress and tangent modulus
    !! follow by automatic differentiation.
    use iso_fortran_env
    use autodiff
    implicit none
    private
    public :: flexible_material
    public :: linear_elastic_1d
    public :: polynomial_elastic_1d

    type, abstract :: flexible_material
        !! Base type for one-dimensional elastic materials.
    contains
        procedure(energy_density_function), deferred :: strain_energy_density
            !! Evaluates the strain energy density for a Green-Lagrange strain.
    end type

    type, extends(flexible_material) :: linear_elastic_1d
        !! Linear elastic material with \(W(E) = \frac{1}{2} M E^2\).
        real(real64) :: modulus = 1.0d0
            !! The elastic modulus.
    contains
        procedure :: strain_energy_density => linear_energy_density
    end type

    type, extends(flexible_material) :: polynomial_elastic_1d
        !! Nonlinear elastic material with
        !! \(W(E) = \frac{1}{2} M E^2 + \frac{1}{3} a E^3 + \frac{1}{4} b E^4\).
        real(real64) :: modulus = 1.0d0
            !! The elastic modulus.
        real(real64) :: quadratic = 0.0d0
            !! The cubic-stress coefficient, \(a\).
        real(real64) :: quartic = 0.0d0
            !! The quartic-stress coefficient, \(b\).
    contains
        procedure :: strain_energy_density => polynomial_energy_density
    end type

    abstract interface
        function energy_density_function(this, e) result(w)
            !! Evaluates the strain energy density of the material for a
            !! Green-Lagrange strain.  The function must be written for
            !! hyper-dual arguments so that the element stress and tangent
            !! follow by automatic differentiation.
            import
            class(flexible_material), intent(in) :: this
                !! The material.
            type(hyperdual), intent(in) :: e
                !! The Green-Lagrange strain.
            type(hyperdual) :: w
                !! The strain energy density.
        end function
    end interface

contains
! ------------------------------------------------------------------------------
function linear_energy_density(this, e) result(w)
    !! Evaluates the strain energy density \(W(E) = \frac{1}{2} M E^2\).
    class(linear_elastic_1d), intent(in) :: this
        !! The material.
    type(hyperdual), intent(in) :: e
        !! The Green-Lagrange strain.
    type(hyperdual) :: w
        !! The strain energy density.

    w = 0.5d0 * this%modulus * e * e
end function

! ------------------------------------------------------------------------------
function polynomial_energy_density(this, e) result(w)
    !! Evaluates the strain energy density
    !! \(W(E) = \frac{1}{2} M E^2 + \frac{1}{3} a E^3 + \frac{1}{4} b E^4\).
    class(polynomial_elastic_1d), intent(in) :: this
        !! The material.
    type(hyperdual), intent(in) :: e
        !! The Green-Lagrange strain.
    type(hyperdual) :: w
        !! The strain energy density.

    w = 0.5d0 * this%modulus * e * e + &
        this%quadratic * e * e * e / 3.0d0 + &
        this%quartic * e * e * e * e / 4.0d0
end function

! ------------------------------------------------------------------------------
end module
