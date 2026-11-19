module dynamics_flexible_bodies_tests
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use autodiff
    use dynamics_flexible_beam_elements
    use dynamics_flexible_body_types, only : fe_node
    use dynamics_flexible_materials
    use dynamics_helper, only : dot_product, matmul
    implicit none
    private
    public :: test_flexible_material
    public :: test_flexible_beam
    public :: test_ad_products
    public :: invalid_ad_product
    public :: invalid_flexible_input

    type, extends(linear_elastic_isotropic_material) :: nonlinear_material
    contains
        procedure :: plane_stress_energy => nonlinear_plane_energy
        procedure :: bending_energy => nonlinear_bending_energy
        procedure :: torsional_energy => nonlinear_torsion_energy
    end type

contains
    pure function nonlinear_plane_energy(this, strain) result(w)
        class(nonlinear_material), intent(in) :: this
        type(hyperdual), intent(in) :: strain
        type(hyperdual) :: w
        w = 0.5d0*this%modulus*strain**2 + this%modulus*strain**4
    end function

    pure function nonlinear_bending_energy(this, ky, kz, iyy, izz) result(w)
        class(nonlinear_material), intent(in) :: this
        type(hyperdual), intent(in) :: ky
        type(hyperdual), intent(in) :: kz
        real(real64), intent(in) :: iyy
        real(real64), intent(in) :: izz
        type(hyperdual) :: w
        w = this%linear_elastic_isotropic_material%bending_energy(ky, kz, iyy, izz) + &
            this%modulus*(iyy*ky**4 + izz*kz**4)
    end function

    pure function nonlinear_torsion_energy(this, twist, j) result(w)
        class(nonlinear_material), intent(in) :: this
        type(hyperdual), intent(in) :: twist
        real(real64), intent(in) :: j
        type(hyperdual) :: w
        w = this%modulus*j*twist**4
    end function

    function test_flexible_material() result(ok)
        logical :: ok
        type(linear_elastic_isotropic_material) :: material
        real(real64), dimension(6,6) :: d
        real(real64), dimension(6) :: strain
        type(dual), dimension(6) :: ds
        type(hyperdual), dimension(6) :: hs
        type(dual) :: de
        type(hyperdual) :: he
        real(real64) :: energy
        real(real64) :: mu
        real(real64) :: lambda
        integer :: i
        material = linear_elastic_isotropic_material(210.0d0, 0.25d0, 2.0d0)
        mu = 84.0d0
        lambda = 84.0d0
        call material%constitutive_matrix(d)
        ok = abs(d(1,1) - (lambda+2*mu)) < 1.0d-12 .and. abs(d(1,2)-lambda) < 1.0d-12
        ok = ok .and. abs(d(4,4)-mu) < 1.0d-12 .and. abs(d(4,5)) < 1.0d-12
        strain = [0.1d0, -0.02d0, 0.03d0, 0.04d0, 0.0d0, 0.0d0]
        do i = 1, 6
            ds(i) = dual(strain(i), merge(1.0d0, 0.0d0, i == 1))
            hs(i) = hyperdual(strain(i), merge(1.0d0, 0.0d0, i == 1), &
                merge(1.0d0, 0.0d0, i == 2), 0.0d0)
        end do
        energy = material%energy_density(strain)
        de = material%energy_density(ds)
        he = material%energy_density(hs)
        ok = ok .and. abs(energy - 0.5d0*dot_product(strain, matmul(d,strain))) < 1.0d-12
        ok = ok .and. abs(de%x-energy) < 1.0d-12 .and. abs(de%dx-dot_product(d(1,:),strain)) < 1.0d-12
        ok = ok .and. abs(he%x-energy) < 1.0d-12 .and. abs(he%dxy-lambda) < 1.0d-12
        he = material%plane_stress_energy(hyperdual(0.1d0, 1.0d0, 1.0d0, 0.0d0))
        ok = ok .and. abs(he%x-1.05d0) < 1.0d-12 .and. abs(he%dx-21.0d0) < 1.0d-12
        ok = ok .and. abs(he%dxy-210.0d0) < 1.0d-12
        he = material%torsional_energy(hyperdual(0.2d0, 1.0d0, 1.0d0, 0.0d0), 0.07d0)
        ok = ok .and. abs(he%dxy-mu*0.07d0) < 1.0d-12
        if (.not.ok) print "(A)", "TEST FAILED: flexible material"
    end function

    function test_ad_products() result(ok)
        logical :: ok
        real(real64), dimension(2,3) :: a
        real(real64), dimension(3) :: x
        real(real64), dimension(3) :: dx
        real(real64), dimension(3) :: dy
        real(real64), dimension(3) :: dxy
        type(dual), dimension(3) :: du
        type(dual), dimension(3) :: dv
        type(dual), dimension(0) :: empty_dual
        type(dual), allocatable, dimension(:) :: dm
        type(dual) :: dd
        type(hyperdual), dimension(3) :: hu
        type(hyperdual), dimension(3) :: hv
        type(hyperdual), dimension(0) :: empty_hyperdual
        type(hyperdual), allocatable, dimension(:) :: hm
        type(hyperdual) :: hd
        integer :: i
        a = reshape([1.0d0, -2.0d0, 3.0d0, 4.0d0, -1.0d0, 2.0d0], [2,3])
        x = [2.0d0, -1.0d0, 3.0d0]
        dx = [1.0d0, 2.0d0, -1.0d0]
        dy = [-2.0d0, 1.0d0, 4.0d0]
        dxy = [3.0d0, -2.0d0, 1.0d0]
        do i = 1, 3
            du(i) = dual(x(i), dx(i))
            dv(i) = dual(dy(i), dxy(i))
            hu(i) = hyperdual(x(i), dx(i), dy(i), dxy(i))
            hv(i) = hyperdual(dy(i), dxy(i), dx(i), x(i))
        end do
        dd = dot_product(du, dv)
        hd = dot_product(hu, hv)
        allocate(dm(2), hm(2))
        dm = matmul(a, du)
        hm = matmul(a, hu)
        ok = abs(dd%x-dot_product(x,dy)) < 1.0d-12
        ok = ok .and. abs(dd%dx-dot_product(dx,dy)-dot_product(x,dxy)) < 1.0d-12
        ok = ok .and. abs(hd%x-dot_product(x,dy)) < 1.0d-12
        ok = ok .and. abs(hd%dx-dot_product(dx,dy)-dot_product(x,dxy)) < 1.0d-12
        ok = ok .and. abs(hd%dy-dot_product(dy,dy)-dot_product(x,dx)) < 1.0d-12
        ok = ok .and. abs(hd%dxy-dot_product(dxy,dy)-dot_product(dx,dx) &
            -dot_product(dy,dxy)-dot_product(x,x)) < 1.0d-12
        ok = ok .and. all(abs(dm%x-matmul(a,x)) < 1.0d-12)
        ok = ok .and. all(abs(dm%dx-matmul(a,dx)) < 1.0d-12)
        ok = ok .and. all(abs(hm%x-matmul(a,x)) < 1.0d-12)
        ok = ok .and. all(abs(hm%dx-matmul(a,dx)) < 1.0d-12)
        ok = ok .and. all(abs(hm%dy-matmul(a,dy)) < 1.0d-12)
        ok = ok .and. all(abs(hm%dxy-matmul(a,dxy)) < 1.0d-12)
        dd = dot_product(empty_dual, empty_dual)
        hd = dot_product(empty_hyperdual, empty_hyperdual)
        ok = ok .and. dd%x == 0.0d0 .and. dd%dx == 0.0d0
        ok = ok .and. hd%x == 0.0d0 .and. hd%dx == 0.0d0 .and. hd%dy == 0.0d0 .and. hd%dxy == 0.0d0
        if (.not.ok) print "(A)", "TEST FAILED: dual and hyperdual products"
    end function

    subroutine invalid_ad_product(selection)
        character(*), intent(in) :: selection
        type(dual), dimension(2) :: d
        type(hyperdual), dimension(2) :: h
        real(real64), dimension(2,3) :: a
        type(dual) :: dd
        type(hyperdual) :: hd
        type(dual), allocatable, dimension(:) :: dm
        type(hyperdual), allocatable, dimension(:) :: hm
        d = dual(1.0d0, 0.0d0)
        h = hyperdual(1.0d0, 0.0d0, 0.0d0, 0.0d0)
        a = 1.0d0
        select case (selection)
        case ("ad-dual-dot")
            dd = dot_product(d, d(1:1))
        case ("ad-hyperdual-dot")
            hd = dot_product(h, h(1:1))
        case ("ad-dual-matmul")
            dm = matmul(a, d)
        case ("ad-hyperdual-matmul")
            hm = matmul(a, h)
        case default
            error stop "Unknown AD product test"
        end select
    end subroutine

    subroutine invalid_flexible_input(selection)
        character(*), intent(in) :: selection
        type(linear_elastic_isotropic_material), target :: material
        type(euler_bernoulli_beam_element) :: beam
        type(fe_node), dimension(2) :: nodes
        real(real64), dimension(12,12) :: m
        real(real64), dimension(6,6) :: d
        real(real64) :: energy
        material%modulus = 100.0d0
        material%poissons_ratio = 0.25d0
        material%density = 1.0d0
        nodes(1) = fe_node(10, [0.0d0, 0.0d0, 0.0d0])
        nodes(2) = fe_node(20, [1.0d0, 0.0d0, 0.0d0])
        beam%material => material
        beam%area = 2.0d0
        beam%j = 0.2d0
        beam%iyy = 0.3d0
        beam%izz = 0.5d0
        beam%node_indices = [1, 2]
        select case (selection)
        case ("flex-invalid-strain")
            energy = material%energy_density([0.0d0])
            print *, energy
        case ("flex-invalid-material")
            material%modulus = ieee_value(0.0d0, ieee_quiet_nan)
            call material%constitutive_matrix(d)
        case ("flex-invalid-section")
            beam%izz = -0.5d0
            call beam%mass_matrix(m, nodes)
        case ("flex-invalid-node-index")
            beam%node_indices = [1, 3]
            call beam%mass_matrix(m, nodes)
        case ("flex-invalid-repeated-node")
            beam%node_indices = [2, 2]
            call beam%mass_matrix(m, nodes)
        case default
            error stop "Unknown flexible input test"
        end select
    end subroutine

    function test_flexible_beam() result(ok)
        logical :: ok
        type(nonlinear_material), target :: material
        type(linear_elastic_isotropic_material), target :: linear_material
        type(euler_bernoulli_beam_element) :: beam
        type(euler_bernoulli_beam_element) :: neighbor
        type(fe_node), dimension(3) :: nodes
        real(real64), dimension(12) :: u
        real(real64), dimension(12) :: f
        real(real64), dimension(12) :: fp
        real(real64), dimension(12) :: fm
        real(real64), dimension(12) :: direction
        real(real64), dimension(12,12) :: k
        real(real64), dimension(12,12) :: m
        type(hyperdual), dimension(12) :: seeded
        type(hyperdual) :: energy
        real(real64) :: strain
        real(real64) :: expected
        real(real64) :: h
        real(real64) :: energy_plus
        real(real64) :: energy_minus
        real(real64) :: s
        real(real64) :: twist
        real(real64), parameter, dimension(4) :: gp = [0.06943184420297371d0, 0.33000947820757187d0, &
            0.66999052179242813d0, 0.93056815579702629d0]
        real(real64), parameter, dimension(4) :: gw = [0.17392742256872693d0, 0.32607257743127307d0, &
            0.32607257743127307d0, 0.17392742256872693d0]
        integer :: i
        material%modulus = 100.0d0
        material%poissons_ratio = 0.25d0
        material%density = 3.0d0
        beam%material => material
        nodes(1) = fe_node(101, [0.0d0, 0.0d0, 0.0d0])
        nodes(2) = fe_node(307, [2.0d0, 0.0d0, 0.0d0])
        nodes(3) = fe_node(502, [4.0d0, 0.0d0, 0.0d0])
        beam%node_indices = [1, 2]
        beam%area = 0.08d0
        beam%j = 0.0007d0
        beam%iyy = 0.0004d0
        beam%izz = 0.0013d0
        beam%equation_indices = [(i, i=1,12)]
        ok = beam%dof_count() == 12 .and. beam%node_count() == 2
        ok = ok .and. beam%dof_per_node() == 6 .and. beam%dimensionality() == 3
        ok = ok .and. all(beam%dof_indices() == [(i, i=1,12)])
        u = 0.0d0
        call beam%internal_force_vector(u, f, nodes)
        call beam%stiffness_matrix(u, k, nodes)
        ok = ok .and. maxval(abs(f)) < 1.0d-11
        ok = ok .and. abs(k(1,1) - material%modulus*beam%area/2.0d0) < 1.0d-10
        ok = ok .and. abs(k(2,2) - 12.0d0*material%modulus*beam%izz/2.0d0**3) < 1.0d-10
        ok = ok .and. abs(k(3,3) - 12.0d0*material%modulus*beam%iyy/2.0d0**3) < 1.0d-10
        linear_material = material%linear_elastic_isotropic_material
        beam%material => linear_material
        call beam%stiffness_matrix(u, k, nodes)
        expected = 0.5d0*material%modulus/(1.0d0+material%poissons_ratio)*beam%j/2.0d0
        ok = ok .and. abs(k(4,4)-expected) < 1.0d-12
        beam%material => material
        call beam%mass_matrix(m, nodes)
        direction = 0.0d0
        direction(1) = 1.0d0
        direction(7) = 1.0d0
        ok = ok .and. abs(dot_product(direction,matmul(m,direction))-material%density*beam%area*2.0d0) < 1.0d-12
        direction = 0.0d0
        direction(4) = 1.0d0
        direction(10) = 1.0d0
        expected = material%density*2.0d0*(beam%iyy+beam%izz)
        ok = ok .and. abs(dot_product(direction,matmul(m,direction))-expected) < 1.0d-12
        ok = ok .and. maxval(abs(m-transpose(m))) < 1.0d-12
        u(7) = 0.4d0
        call beam%internal_force_vector(u, f, nodes)
        strain = 0.5d0*(1.2d0**2-1.0d0)
        expected = beam%area*material%modulus*(strain+4.0d0*strain**3)*1.2d0
        ok = ok .and. abs(f(7)-expected) < 1.0d-10 .and. abs(f(1)+expected) < 1.0d-10
        u = 0.0d0
        u(7) = 2.0d0*(cos(0.8d0)-1.0d0)
        u(8) = 2.0d0*sin(0.8d0)
        u(6) = 0.8d0
        u(12) = 0.8d0
        call beam%internal_force_vector(u, f, nodes)
        ok = ok .and. maxval(abs(f)) < 1.0d-10
        u = 0.0d0
        u(10) = 0.3d0
        do i = 1, 12
            seeded(i) = hyperdual(u(i), 0.0d0, 0.0d0, 0.0d0)
        end do
        energy = beam%strain_energy(seeded, nodes)
        expected = 0.0d0
        do i = 1, 4
            s = gp(i)
            twist = sin(0.3d0) / (2.0d0*((1.0d0-s+s*cos(0.3d0))**2+(s*sin(0.3d0))**2))
            expected = expected + gw(i)*2.0d0*material%modulus*beam%j*twist**4
        end do
        ok = ok .and. abs(energy%x-expected) < 1.0d-12
        u = [0.01d0, 0.03d0, -0.02d0, 0.04d0, -0.03d0, 0.05d0, &
             0.12d0, -0.04d0, 0.06d0, -0.02d0, 0.03d0, -0.01d0]
        direction = [0.3d0, -0.2d0, 0.1d0, 0.4d0, -0.3d0, 0.2d0, &
                     -0.1d0, 0.2d0, 0.3d0, -0.4d0, 0.1d0, 0.2d0]
        h = 1.0d-5
        call beam%internal_force_vector(u, f, nodes)
        call beam%stiffness_matrix(u, k, nodes)
        call beam%internal_force_vector(u+h*direction, fp, nodes)
        call beam%internal_force_vector(u-h*direction, fm, nodes)
        do i = 1, 12
            seeded(i) = hyperdual(u(i)+h*direction(i), 0.0d0, 0.0d0, 0.0d0)
        end do
        energy = beam%strain_energy(seeded, nodes)
        energy_plus = energy%x
        do i = 1, 12
            seeded(i) = hyperdual(u(i)-h*direction(i), 0.0d0, 0.0d0, 0.0d0)
        end do
        energy = beam%strain_energy(seeded, nodes)
        energy_minus = energy%x
        ok = ok .and. abs((energy_plus-energy_minus)/(2*h)-dot_product(f,direction)) < 1.0d-8
        ok = ok .and. maxval(abs(k-transpose(k))) < 1.0d-10
        ok = ok .and. maxval(abs((fp-fm)/(2*h)-matmul(k,direction))) < 1.0d-7

        ! Both elements read node 2 from the same table, regardless of its ID.
        neighbor = beam
        neighbor%node_indices = [2, 3]
        u = 0.0d0
        call neighbor%stiffness_matrix(u, k, nodes)
        ok = ok .and. abs(k(1,1) - material%modulus*beam%area/2.0d0) < 1.0d-10
        nodes(2)%coordinates(1) = 3.0d0
        call beam%stiffness_matrix(u, k, nodes)
        ok = ok .and. abs(k(1,1) - material%modulus*beam%area/3.0d0) < 1.0d-10
        call neighbor%stiffness_matrix(u, k, nodes)
        ok = ok .and. abs(k(1,1) - material%modulus*beam%area) < 1.0d-10
        if (.not.ok) print "(A)", "TEST FAILED: flexible beam"
    end function
end module
