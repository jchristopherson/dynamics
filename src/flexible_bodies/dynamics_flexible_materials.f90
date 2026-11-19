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
    !! Material models for flexible-body elements. A material supplies an
    !! axial energy law and the moduli for quadratic bending, torsion and shear.
    use iso_fortran_env, only : real64
    use ieee_arithmetic, only : ieee_is_finite
    use autodiff
    implicit none
    private
    public :: flexible_material
    public :: linear_elastic
    public :: polynomial_elastic

    type, abstract :: flexible_material
        !! Base type for flexible-body elastic materials.
        !! Trusses use the axial energy law directly. The default beam section
        !! law also uses the elastic moduli to convert section geometry into
        !! bending, torsion, and shear stiffnesses. Section geometry is not
        !! part of the material model.
        real(real64) :: density = 1.0d0
            !! Reference mass per unit volume. Mass integration requires a finite, nonnegative value.
    contains
        procedure(energy_density_function), deferred :: strain_energy_density
            !! Evaluates stored energy per unit reference volume for the supplied axial strain measure.
        procedure(elastic_moduli_subroutine), deferred :: get_elastic_moduli
            !! Returns the strain-independent Young's and shear moduli used by the default beam section law.
    end type

    type, extends(flexible_material) :: linear_elastic
        !! Isotropic linear material with axial energy per unit reference
        !! volume \(W(e) = \frac{1}{2} E e^2\). The shear modulus is derived
        !! from \(E\) and \(\nu\); it does not enter the axial energy.
        real(real64) :: youngs_modulus = 1.0d0
            !! Young's modulus \(E\), in stress units.
        real(real64) :: poissons_ratio = 0.3d0
            !! Dimensionless Poisson's ratio \(\nu\); the moduli routine requires \(-1 < \nu < 0.5\).
    contains
        procedure :: strain_energy_density => linear_energy_density
        procedure :: get_elastic_moduli => isotropic_elastic_moduli
    end type

    type, extends(linear_elastic) :: polynomial_elastic
        !! Nonlinear axial energy per unit reference volume
        !! \(W(e) = \frac{1}{2} E e^2 + \frac{1}{3} a e^3 + \frac{1}{4} b e^4\).
        !! The inherited \(E\) is the axial tangent at zero strain. The
        !! default beam section law still uses constant \(E\) and \(G\) for
        !! bending, torsion, and shear. Setting \(a=b=0\) recovers the linear law.
        real(real64) :: quadratic = 0.0d0
            !! Coefficient \(a\) of the quadratic term in \(dW/de\), in stress units.
        real(real64) :: quartic = 0.0d0
            !! Coefficient \(b\) of the cubic term in \(dW/de\), in stress units.
    contains
        procedure :: strain_energy_density => polynomial_energy_density
    end type

    interface linear_elastic
        !! Constructs an isotropic linear material from optional elastic constants and density.
        module procedure :: make_linear_elastic
    end interface

    interface polynomial_elastic
        !! Constructs a polynomial material from optional elastic constants, coefficients, and density.
        module procedure :: make_polynomial_elastic
    end interface

    abstract interface
        function energy_density_function(this, e) result(w)
            !! Evaluates stored energy for a scalar axial strain measure.
            !! Implementations must preserve the hyper-dual derivative components
            !! so the element can differentiate the energy for forces and tangent stiffness.
            import
            class(flexible_material), intent(in) :: this
                !! Material model providing the axial constitutive law.
            type(hyperdual), intent(in) :: e
                !! Dimensionless axial strain supplied by the caller: Green-Lagrange
                !! strain for trusses or generalized extension for the default beam section law.
            type(hyperdual) :: w
                !! Stored axial strain energy per unit reference volume.
        end function

        subroutine elastic_moduli_subroutine(this, e, g)
            !! Supplies constant elastic moduli for quadratic section terms.
            !! The default beam section law uses these independently of the current
            !! axial strain and its tangent.
            import :: flexible_material, real64
            class(flexible_material), intent(in) :: this
                !! Material model providing the elastic constants.
            real(real64), intent(out) :: e
                !! Finite positive Young's modulus \(E\), in stress units.
            real(real64), intent(out) :: g
                !! Finite positive shear modulus \(G\), in stress units.
        end subroutine
    end interface

contains
function make_linear_elastic(youngs_modulus, poissons_ratio, density) result(material)
    !! Keeps the original positional order of the elastic constants. Density
    !! is optional and follows them, although inherited storage comes first.
    real(real64), intent(in), optional :: youngs_modulus
        !! Young's modulus in stress units; defaults to 1.
    real(real64), intent(in), optional :: poissons_ratio
        !! Dimensionless Poisson's ratio; defaults to 0.3.
    real(real64), intent(in), optional :: density
        !! Reference mass per unit volume; defaults to 1.
    type(linear_elastic) :: material
        !! Material initialized with the supplied values and defaults.

    if (present(youngs_modulus)) material%youngs_modulus = youngs_modulus
    if (present(poissons_ratio)) material%poissons_ratio = poissons_ratio
    if (present(density)) material%density = density
end function

function make_polynomial_elastic(youngs_modulus, poissons_ratio, quadratic, quartic, density) result(material)
    !! Preserves the flattened constructor's original positional order. A
    !! linear_elastic parent can still be passed to the structure constructor.
    real(real64), intent(in), optional :: youngs_modulus
        !! Axial tangent at zero strain, in stress units; defaults to 1.
    real(real64), intent(in), optional :: poissons_ratio
        !! Dimensionless Poisson's ratio; defaults to 0.3.
    real(real64), intent(in), optional :: quadratic
        !! Coefficient of the quadratic term in the axial stress; defaults to zero.
    real(real64), intent(in), optional :: quartic
        !! Coefficient of the cubic term in the axial stress; defaults to zero.
    real(real64), intent(in), optional :: density
        !! Reference mass per unit volume; defaults to 1.
    type(polynomial_elastic) :: material
        !! Material initialized with the supplied values and defaults.

    material%linear_elastic = make_linear_elastic(youngs_modulus, poissons_ratio, density)
    if (present(quadratic)) material%quadratic = quadratic
    if (present(quartic)) material%quartic = quartic
end function

! ------------------------------------------------------------------------------
subroutine isotropic_elastic_moduli(this, e, g)
    !! Returns \(E\) and computes \(G = E/[2(1+\nu)]\) from isotropic elasticity.
    !! Stops with an error if \(\nu\) is nonfinite or outside the open interval
    !! \((-1, 0.5)\). The default beam section law checks the returned moduli.
    class(linear_elastic), intent(in) :: this
        !! Isotropic material holding Young's modulus and Poisson's ratio.
    real(real64), intent(out) :: e
        !! Finite positive Young's modulus \(E\), in stress units.
    real(real64), intent(out) :: g
        !! Finite positive shear modulus \(G\), in stress units.
    real(real64) :: nu

    e = this%youngs_modulus
    nu = this%poissons_ratio
    if (.not. ieee_is_finite(nu)) error stop "Poisson's ratio must be finite"
    if (nu <= -1.0d0 .or. nu >= 0.5d0) error stop "Poisson's ratio must be between -1 and 0.5"
    g = e / (2.0d0 * (1.0d0 + nu))
end subroutine

! ------------------------------------------------------------------------------
function linear_energy_density(this, e) result(w)
    !! Evaluates the quadratic axial energy \(W(e) = \frac{1}{2} E e^2\).
    !! Poisson's ratio does not enter the axial law.
    class(linear_elastic), intent(in) :: this
        !! Material supplying Young's modulus \(E\).
    type(hyperdual), intent(in) :: e
        !! Dimensionless hyper-dual axial strain in the caller's strain measure.
    type(hyperdual) :: w
        !! Stored axial strain energy per unit reference volume.

    w = 0.5d0 * this%youngs_modulus * e * e
end function

! ------------------------------------------------------------------------------
function polynomial_energy_density(this, e) result(w)
    !! Evaluates \(W(e) = \frac{1}{2} E e^2 + \frac{1}{3} a e^3 + \frac{1}{4} b e^4\).
    !! The conjugate axial stress is \(dW/de = E e + a e^2 + b e^3\).
    class(polynomial_elastic), intent(in) :: this
        !! Material supplying \(E\) and the nonlinear stress coefficients \(a\) and \(b\).
    type(hyperdual), intent(in) :: e
        !! Dimensionless hyper-dual axial strain in the caller's strain measure.
    type(hyperdual) :: w
        !! Stored axial strain energy per unit reference volume.

    w = 0.5d0 * this%youngs_modulus * e * e + &
        this%quadratic * e * e * e / 3.0d0 + &
        this%quartic * e * e * e * e / 4.0d0
end function

! ------------------------------------------------------------------------------
end module
