module dynamics_flexible_materials
    !! Material energy interfaces for flexible elements.
    !! Use a consistent unit system: energy density is energy per reference
    !! volume, torsional section energy is energy per reference length,
    !! moduli are stress, and density is mass per reference volume.
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use autodiff
    use dynamics_helper
    implicit none
    private
    public :: flexible_material
    public :: linear_elastic_isotropic_material

    type, abstract :: flexible_material
        !! Constitutive energy contract for flexible elements. Implementations
        !! provide consistent real, dual and hyper-dual 3D energy evaluations
        !! and separate uniaxial plane-stress and rectangular-section torsion
        !! laws for beams.
        real(real64) :: density
            !! Material density per reference volume.
    contains
        procedure(get_real_energy_density), deferred, public, pass :: real_energy_density
            !! Real-valued 3D strain energy per reference volume.
        procedure(get_dual_energy_density), deferred, public, pass :: dual_energy_density
            !! The same energy with first-order dual derivatives.
        procedure(get_hyperdual_energy_density), deferred, public, pass :: hyperdual_energy_density
            !! The same energy with hyper-dual derivatives.
        procedure(get_plane_stress_energy), deferred, public, pass :: plane_stress_energy
            !! Uniaxial fiber energy with free transverse stresses.
        procedure(get_torsional_energy), deferred, public, pass :: torsional_energy
            !! Rectangular-section torsional energy per reference length.
        generic, public :: energy_density => real_energy_density, &
            dual_energy_density, hyperdual_energy_density
            !! Dispatch the 3D energy evaluation by scalar type.
    end type
    
    abstract interface
        pure function get_plane_stress_energy(this, strain) result(rst)
            !! Evaluate uniaxial fiber energy per reference volume after
            !! imposing the material's free-transverse-stress response.
            !! Derivatives must propagate through the hyper-dual input.
            import flexible_material, hyperdual
            class(flexible_material), intent(in) :: this
                !! Material supplying the reduced constitutive law.
            type(hyperdual), intent(in) :: strain
                !! Dimensionless axial Green-Lagrange strain and seeds.
            type(hyperdual) :: rst
                !! Energy per reference volume and seeded derivatives.
        end function

        pure function get_torsional_energy(this, twist, width, height) result(rst)
            !! Evaluate the rectangular section's torsional energy per
            !! reference length. The material implementation owns the torsion
            !! law and must propagate all hyper-dual derivatives of twist.
            import flexible_material, hyperdual, real64
            class(flexible_material), intent(in) :: this
                !! Material supplying the section torsion law.
            type(hyperdual), intent(in) :: twist
                !! Twist rate per reference length and derivative seeds.
            real(real64), intent(in) :: width
                !! Undeformed section dimension along local y.
            real(real64), intent(in) :: height
                !! Undeformed section dimension along local z.
            type(hyperdual) :: rst
                !! Torsional energy per reference length and derivatives.
        end function

        pure function get_real_energy_density(this, strain) result(rst)
            !! Evaluate 3D strain energy per reference volume.
            use iso_fortran_env, only : real64
            import flexible_material
            class(flexible_material), intent(in) :: this
                !! The [[flexible_material]] type.
            real(real64), intent(in), dimension(:) :: strain
                !! Strains [xx, yy, zz, xy, yz, zx], with engineering shear strains.
            real(real64) :: rst
                !! Energy per reference volume.
        end function

        pure function get_dual_energy_density(this, strain) result(rst)
            !! Evaluate 3D strain energy and first derivatives.
            use autodiff, only : dual
            import flexible_material
            class(flexible_material), intent(in) :: this
                !! The [[flexible_material]] type.
            type(dual), intent(in), dimension(:) :: strain
                !! Strains [xx, yy, zz, xy, yz, zx], with engineering shear strains.
            type(dual) :: rst
                !! Energy per reference volume and its directional derivative.
        end function

        pure function get_hyperdual_energy_density(this, strain) result(rst)
            !! Evaluate 3D strain energy and first and mixed derivatives.
            use autodiff, only : hyperdual
            import flexible_material
            class(flexible_material), intent(in) :: this
                !! The [[flexible_material]] type.
            type(hyperdual), intent(in), dimension(:) :: strain
                !! Strains [xx, yy, zz, xy, yz, zx], with engineering shear strains.
            type(hyperdual) :: rst
                !! Energy per reference volume and its seeded derivatives.
        end function
    end interface

    type, extends(flexible_material) :: linear_elastic_isotropic_material
        !! Isotropic quadratic strain energy with Young's modulus and
        !! Poisson's ratio. The beam plane-stress law is E*strain**2/2;
        !! torsion uses the rectangular Saint-Venant torsion constant.
        real(real64) :: modulus
            !! The material modulus of elasticity.
        real(real64) :: poissons_ratio = 0.3d0
            !! The material Poisson's ratio.
    contains
        procedure, public :: constitutive_matrix => lei_constitutive_matrix
        procedure, public :: real_energy_density => lei_real_energy_density
        procedure, public :: dual_energy_density => lei_dual_energy_density
        procedure, public :: hyperdual_energy_density => lei_hyperdual_energy_density
        procedure, public :: plane_stress_energy => lei_plane_stress_energy
        procedure, public :: torsional_energy => lei_torsional_energy
    end type

    interface linear_elastic_isotropic_material
        module procedure :: lei_initialize
    end interface

contains
! ------------------------------------------------------------------------------
    pure function lei_initialize(e, nu, rho) result(rst)
        !! Construct an isotropic material from the supplied properties.
        !! Properties are checked by the energy and constitutive routines,
        !! not by this constructor.
        real(real64), intent(in) :: e
            !! Young's modulus, in units of stress.
        real(real64), intent(in) :: nu
            !! Dimensionless Poisson's ratio.
        real(real64), intent(in) :: rho
            !! Density per reference volume.
        type(linear_elastic_isotropic_material) :: rst
            !! The [[linear_elastic_isotropic_material]] object.
        rst%modulus = e
        rst%poissons_ratio = nu
        rst%density = rho
    end function

! ------------------------------------------------------------------------------
    pure subroutine lei_constitutive_matrix(this, d)
        !! Compute the isotropic 3D constitutive matrix in the order
        !! [xx, yy, zz, xy, yz, zx], using engineering shear strains.
        !! Requires positive Young's modulus and -1 < Poisson's ratio < 0.5;
        !! otherwise execution terminates with an error.
        class(linear_elastic_isotropic_material), intent(in) :: this
            !! The [[linear_elastic_isotropic_material]] object.
        real(real64), intent(out), dimension(6,6) :: d
            !! The 6-by-6 stress-strain matrix.

        real(real64) :: lambda
        real(real64) :: mu
        integer :: i
        integer :: j

        if (.not.ieee_is_finite(this%modulus) .or. .not.ieee_is_finite(this%poissons_ratio) .or. &
            this%modulus <= 0.0d0 .or. this%poissons_ratio <= -1.0d0 .or. &
            this%poissons_ratio >= 0.5d0) error stop "Invalid isotropic material properties"
        mu = this%modulus / (2.0d0 * (1.0d0 + this%poissons_ratio))
        lambda = this%modulus * this%poissons_ratio / &
            ((1.0d0 + this%poissons_ratio) * (1.0d0 - 2.0d0 * this%poissons_ratio))
        d = 0.0d0
        do i = 1, 3
            do j = 1, 3
                d(i,j) = lambda
            end do
            d(i,i) = lambda + 2.0d0 * mu
            d(i+3,i+3) = mu
        end do
    end subroutine

! ------------------------------------------------------------------------------
    pure function lei_plane_stress_energy(this, strain) result(rst)
        !! Evaluate E*strain**2/2 for a uniaxial fiber with free transverse
        !! stresses. A nonpositive Young's modulus terminates with an error.
        class(linear_elastic_isotropic_material), intent(in) :: this
            !! Material supplying Young's modulus.
        type(hyperdual), intent(in) :: strain
            !! Axial Green-Lagrange strain and derivative seeds.
        type(hyperdual) :: rst
            !! Energy per reference volume with propagated derivatives.
        if (.not.ieee_is_finite(this%modulus) .or. this%modulus <= 0.0d0) error stop "Invalid Young's modulus"
        rst = 0.5d0 * this%modulus * strain * strain
    end function

! ------------------------------------------------------------------------------
    pure function lei_torsional_energy(this, twist, width, height) result(rst)
        !! Evaluate G*J*twist**2/2 per reference length, where G is the
        !! isotropic shear modulus and J is the rectangular Saint-Venant
        !! torsion constant. J is approximated by 200 terms of its odd-index
        !! series. Invalid section dimensions or elastic properties terminate
        !! with an error.
        class(linear_elastic_isotropic_material), intent(in) :: this
            !! Material supplying Young's modulus and Poisson's ratio.
        type(hyperdual), intent(in) :: twist
            !! Twist rate per reference length and derivative seeds.
        real(real64), intent(in) :: width
            !! Positive section dimension along local y.
        real(real64), intent(in) :: height
            !! Positive section dimension along local z.
        type(hyperdual) :: rst
            !! Torsional energy per reference length and its derivatives.
        real(real64) :: a
        real(real64) :: b
        real(real64) :: pi
        real(real64) :: series
        real(real64) :: torsion_constant
        real(real64) :: n
        integer :: i

        if (.not.ieee_is_finite(this%modulus) .or. .not.ieee_is_finite(this%poissons_ratio) .or. &
            this%modulus <= 0.0d0 .or. this%poissons_ratio <= -1.0d0 .or. &
            this%poissons_ratio >= 0.5d0) error stop "Invalid isotropic material properties"
        if (.not.ieee_is_finite(width) .or. .not.ieee_is_finite(height) .or. &
            width <= 0.0d0 .or. height <= 0.0d0) error stop "Invalid rectangular section"
        a = max(width, height)
        b = min(width, height)
        pi = acos(-1.0d0)
        series = 0.0d0
        do i = 1, 200
            n = real(2 * i - 1, real64)
            series = series + tanh(n * pi * a / (2.0d0 * b)) / n**5
        end do
        torsion_constant = a * b**3 * (1.0d0 / 3.0d0 - 64.0d0 * b * series / (pi**5 * a))
        rst = 0.5d0 * (this%modulus / (2.0d0 * (1.0d0 + this%poissons_ratio))) * &
            torsion_constant * twist * twist
    end function

! ------------------------------------------------------------------------------
    pure function lei_real_energy_density(this, strain) result(rst)
        !! Compute quadratic isotropic 3D energy per reference volume.
        !! Invalid material properties or a strain vector whose length is not
        !! six terminate with an error.
        class(linear_elastic_isotropic_material), intent(in) :: this
            !! The [[linear_elastic_isotropic_material]] type.
        real(real64), intent(in), dimension(:) :: strain
            !! [xx, yy, zz, xy, yz, zx], using engineering shear strains.
        real(real64) :: rst
            !! The strain energy density.

        real(real64) :: D(6, 6)
        if (size(strain) /= 6) error stop "Material strain must have six components"
        call this%constitutive_matrix(D)
        rst = 0.5d0 * dot_product(strain, matmul(D, strain))
    end function

! ------------------------------------------------------------------------------
    pure function lei_dual_energy_density(this, strain) result(rst)
        !! Compute quadratic isotropic 3D energy and its dual derivative.
        !! Invalid material properties or a strain vector whose length is not
        !! six terminate with an error.
        class(linear_elastic_isotropic_material), intent(in) :: this
            !! The [[linear_elastic_isotropic_material]] type.
        type(dual), intent(in), dimension(:) :: strain
            !! [xx, yy, zz, xy, yz, zx], using engineering shear strains.
        type(dual) :: rst
            !! Energy per reference volume and its directional derivative.

        real(real64) :: D(6, 6)
        if (size(strain) /= 6) error stop "Material strain must have six components"
        call this%constitutive_matrix(D)
        rst = 0.5d0 * dot_product(strain, matmul(D, strain))
    end function

! ------------------------------------------------------------------------------
    pure function lei_hyperdual_energy_density(this, strain) result(rst)
        !! Compute quadratic isotropic 3D energy and its hyper-dual
        !! derivatives. Invalid material properties or a strain vector whose
        !! length is not six terminate with an error.
        class(linear_elastic_isotropic_material), intent(in) :: this
            !! The [[linear_elastic_isotropic_material]] type.
        type(hyperdual), intent(in), dimension(:) :: strain
            !! [xx, yy, zz, xy, yz, zx], using engineering shear strains.
        type(hyperdual) :: rst
            !! Energy per reference volume and its seeded derivatives.

        real(real64) :: D(6, 6)
        if (size(strain) /= 6) error stop "Material strain must have six components"
        call this%constitutive_matrix(D)
        rst = 0.5d0 * dot_product(strain, matmul(D, strain))
    end function

! ------------------------------------------------------------------------------
end module