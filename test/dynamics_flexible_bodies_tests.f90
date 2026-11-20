module dynamics_flexible_bodies_tests
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use autodiff
    use dynamics_flexible_beam_elements
    use dynamics_flexible_body_types, only : fe_node, degree_of_freedom, DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ
    use dynamics_flexible_assembly, only : finite_element_mesh
    use dynamics_boundary_conditions, only : nodal_condition, &
        apply_force_boundary_conditions, apply_displacement_boundary_conditions, &
        restore_displacement_boundary_conditions, resolve_boundary_conditions
    use dynamics_flexible_body_solvers, only : flexible_static_settings, flexible_static_info, solve_flexible_static
    use linalg, only : csr_matrix, csr_to_dense
    use dynamics_flexible_materials
    use dynamics_helper, only : dot_product, matmul
    implicit none
    private
    public :: test_flexible_material
    public :: test_flexible_beam
    public :: test_flexible_mesh
    public :: test_flexible_mesh_boundary_conditions
    public :: test_flexible_static
    public :: test_ad_products
    public :: invalid_ad_product
    public :: invalid_flexible_input

    type, extends(linear_elastic_isotropic_material) :: nonlinear_material
    contains
        procedure :: plane_stress_energy => nonlinear_plane_energy
        procedure :: bending_energy => nonlinear_bending_energy
        procedure :: torsional_energy => nonlinear_torsion_energy
    end type

    ! Topology-only test double; numerical beam routines are not used.
    type, extends(euler_bernoulli_beam_element) :: labeled_test_beam
        integer, dimension(6) :: labels = [DOF_UY, DOF_UX, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ]
    contains
        procedure :: dof_components => test_beam_components
    end type

    type, extends(euler_bernoulli_beam_element) :: static_test_element
    contains
        procedure :: internal_force_vector => static_test_force
        procedure :: stiffness_matrix => static_test_stiffness
    end type

contains
    subroutine static_test_force(this, q, f, nodes)
        class(static_test_element), intent(inout), target :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: f
        type(fe_node), intent(in), dimension(:) :: nodes
        f = q + q**3
        f(:size(q)-1) = f(:size(q)-1) + 0.2d0*q(2:)
        f(2:) = f(2:) + 0.2d0*q(:size(q)-1)
    end subroutine

    subroutine static_test_stiffness(this, q, k, nodes)
        class(static_test_element), intent(inout), target :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: k
        type(fe_node), intent(in), target, dimension(:) :: nodes
        integer :: i
        k = 0.0d0
        do i = 1, size(q)
            k(i,i) = 1.0d0 + 3.0d0*q(i)**2
            if (i < size(q)) then
                k(i,i+1) = 0.2d0
                k(i+1,i) = 0.2d0
            end if
        end do
    end subroutine

    function test_flexible_static() result(ok)
        logical :: ok
        type(nonlinear_material), target :: material
        type(linear_elastic_isotropic_material), target :: linear_material
        type(euler_bernoulli_beam_element) :: beam
        type(static_test_element) :: element
        type(finite_element_mesh) :: mesh
        type(finite_element_mesh) :: large_mesh
        type(flexible_static_settings) :: settings
        type(flexible_static_info) :: info
        type(nodal_condition), dimension(11) :: conditions
        type(nodal_condition), dimension(0) :: no_conditions
        type(nodal_condition), dimension(12) :: all_conditions
        type(nodal_condition), dimension(1) :: one_condition
        real(real64), dimension(12) :: q
        real(real64), dimension(12) :: load
        real(real64), dimension(12) :: before
        real(real64), dimension(48) :: large_q
        real(real64), dimension(48) :: target_q
        real(real64), dimension(48) :: large_load
        real(real64) :: strain
        integer :: i
        integer :: method

        material%modulus = 100.0d0
        material%poissons_ratio = 0.25d0
        material%density = 1.0d0
        beam%material => material
        beam%area = 1.0d0
        beam%j = 0.1d0
        beam%iyy = 0.2d0
        beam%izz = 0.3d0
        beam%node_indices = [1, 2]
        mesh%nodes = [fe_node(1, [0.0d0, 0.0d0, 0.0d0]), fe_node(2, [1.0d0, 0.0d0, 0.0d0])]
        allocate(mesh%elements(1))
        allocate(mesh%elements(1)%element, source=beam)
        conditions = [nodal_condition(1, DOF_UX, 0.1d0), nodal_condition(1, DOF_UY), &
            nodal_condition(1, DOF_UZ), nodal_condition(1, DOF_RX), nodal_condition(1, DOF_RY), &
            nodal_condition(1, DOF_RZ), nodal_condition(2, DOF_UY), nodal_condition(2, DOF_UZ), &
            nodal_condition(2, DOF_RX), nodal_condition(2, DOF_RY), nodal_condition(2, DOF_RZ)]
        strain = (1.3d0**2 - 1.0d0) / 2.0d0
        load = 0.0d0
        load(7) = 100.0d0 * (strain + 4.0d0*strain**3) * 1.3d0
        ok = .true.
        do method = 1, 2
            settings%sparse = method == 2
            q = 0.0d0
            call solve_flexible_static(mesh, load, conditions, q, settings, info)
            ok = ok .and. info%converged .and. info%iterations > 1
            ok = ok .and. abs(q(1)-0.1d0) < 1.0d-12 .and. abs(q(7)-0.4d0) < 1.0d-8
            ok = ok .and. info%residual_norm <= settings%absolute_tolerance + &
                settings%relative_tolerance * abs(load(7))
        end do

        ! The same finite-strain equilibrium with a linear material.
        linear_material = linear_elastic_isotropic_material(100.0d0, 0.25d0, 1.0d0)
        select type (mesh_beam => mesh%elements(1)%element)
        type is (euler_bernoulli_beam_element)
            mesh_beam%material => linear_material
        end select
        load(7) = 100.0d0 * strain * 1.3d0
        settings%sparse = .false.
        q = 0.0d0
        call solve_flexible_static(mesh, load, conditions, q, settings, info)
        ok = ok .and. info%converged .and. abs(q(7)-0.4d0) < 1.0d-8
        select type (mesh_beam => mesh%elements(1)%element)
        type is (euler_bernoulli_beam_element)
            mesh_beam%material => material
        end select
        load(7) = 100.0d0 * (strain + 4.0d0*strain**3) * 1.3d0

        ! A failed solve must not publish the trial state.
        settings%max_iterations = 1
        q = 0.0d0
        before = q
        call solve_flexible_static(mesh, load, conditions, q, settings, info)
        ok = ok .and. .not.info%converged .and. all(q == before)
        settings%max_iterations = 30
        settings%sparse = .true.
        q = 0.0d0
        call solve_flexible_static(mesh, load, no_conditions, q, settings, info)
        ok = ok .and. .not.info%converged .and. info%linear_solve_failed .and. all(q == 0.0d0)
        all_conditions(:11) = conditions
        all_conditions(12) = nodal_condition(2, DOF_UX, 0.4d0)
        call solve_flexible_static(mesh, load, all_conditions, q, settings, info)
        ok = ok .and. info%converged .and. abs(q(7)-0.4d0) < 1.0d-12

        ! Exercise the actual sparse iterative path above 32 free DOFs.
        large_mesh = finite_element_mesh(8, 4)
        do i = 1, 8
            large_mesh%nodes(i) = fe_node(i, [real(i,real64), 0.0d0, 0.0d0])
        end do
        do i = 1, 4
            element%node_indices = [2*i-1, 2*i]
            call large_mesh%assign_element(i, element)
        end do
        target_q = [(0.01d0*i, i=1,48)]
        large_load = target_q + target_q**3
        do i = 1, 4
            large_load(12*i-11:12*i-1) = large_load(12*i-11:12*i-1) + 0.2d0*target_q(12*i-10:12*i)
            large_load(12*i-10:12*i) = large_load(12*i-10:12*i) + 0.2d0*target_q(12*i-11:12*i-1)
        end do
        large_q = 0.0d0
        call solve_flexible_static(large_mesh, large_load, no_conditions, large_q, settings, info)
        ok = ok .and. info%converged .and. maxval(abs(large_q-target_q)) < 1.0d-8

        ! Full Newton overshoots here; backtracking must recover.
        target_q = 2.0d0
        large_load = target_q + target_q**3
        do i = 1, 4
            large_load(12*i-11:12*i-1) = large_load(12*i-11:12*i-1) + 0.2d0*target_q(12*i-10:12*i)
            large_load(12*i-10:12*i) = large_load(12*i-10:12*i) + 0.2d0*target_q(12*i-11:12*i-1)
        end do
        settings%max_backtracks = 0
        large_q = 0.0d0
        call solve_flexible_static(large_mesh, large_load, no_conditions, large_q, settings, info)
        ok = ok .and. .not.info%converged .and. all(large_q == 0.0d0)
        settings%max_backtracks = 20
        call solve_flexible_static(large_mesh, large_load, no_conditions, large_q, settings, info)
        ok = ok .and. info%converged .and. maxval(abs(large_q-target_q)) < 1.0d-7
        settings%sparse = .false.
        large_q = 0.0d0
        call solve_flexible_static(large_mesh, large_load, no_conditions, large_q, settings, info)
        ok = ok .and. info%converged .and. maxval(abs(large_q-target_q)) < 1.0d-7
        settings%sparse = .true.
        one_condition(1) = nodal_condition(1, DOF_UX, 2.0d0)
        large_q = 0.0d0
        call solve_flexible_static(large_mesh, large_load, one_condition, large_q, settings, info)
        ok = ok .and. info%converged .and. maxval(abs(large_q-target_q)) < 1.0d-7
        call solve_flexible_static(large_mesh, large_load, one_condition, large_q, settings, info)
        ok = ok .and. info%converged .and. info%iterations == 0
        if (.not.ok) print "(A)", "TEST FAILED: nonlinear flexible static solver"
    end function

    pure function test_beam_components(this) result(components)
        class(labeled_test_beam), intent(in) :: this
        integer, allocatable, dimension(:) :: components

        components = this%labels
    end function

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
        type(labeled_test_beam) :: labeled
        type(finite_element_mesh) :: mesh
        type(flexible_static_settings) :: settings
        integer, allocatable, dimension(:) :: offsets
        integer, allocatable, dimension(:) :: equations
        type(fe_node), dimension(2) :: nodes
        real(real64), dimension(12,12) :: m
        real(real64), dimension(6,6) :: d
        real(real64), dimension(12) :: q
        real(real64), dimension(12) :: load
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
        q = 0.0d0
        load = 0.0d0
        select case (selection)
        case ("flex-invalid-static-duplicate")
            mesh%nodes = nodes
            allocate(mesh%elements(1))
            allocate(mesh%elements(1)%element, source=beam)
            call solve_flexible_static(mesh, load, &
                [nodal_condition(1, DOF_UX), nodal_condition(1, DOF_UX)], q)
        case ("flex-invalid-static-options")
            settings%relative_tolerance = -1.0d0
            call solve_flexible_static(mesh, load, [nodal_condition::], q, settings)
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
        case ("flex-invalid-mesh-node")
            mesh%nodes = nodes
            beam%node_indices = [1, 3]
            allocate(mesh%elements(1))
            allocate(mesh%elements(1)%element, source=beam)
            call mesh%number_dofs(offsets)
        case ("flex-invalid-dof-order")
            mesh%nodes = [nodes, fe_node(30, [2.0d0, 0.0d0, 0.0d0])]
            allocate(mesh%elements(2))
            allocate(mesh%elements(1)%element, source=beam)
            labeled%node_indices = [2, 3]
            allocate(mesh%elements(2)%element, source=labeled)
            call mesh%number_dofs(offsets)
        case ("flex-invalid-dof-label")
            mesh%nodes = nodes
            allocate(mesh%elements(1))
            labeled%node_indices = [1, 2]
            labeled%labels(2) = 7
            allocate(mesh%elements(1)%element, source=labeled)
            call resolve_boundary_conditions(mesh, [nodal_condition(1, DOF_UX)], equations)
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
   function test_flexible_mesh() result(ok)
       logical :: ok
       type(linear_elastic_isotropic_material), target :: material
       type(euler_bernoulli_beam_element) :: beam
       type(finite_element_mesh) :: mesh
       type(csr_matrix) :: sparse
       type(csr_matrix) :: threaded_sparse
       integer, allocatable, dimension(:) :: offsets
       integer, dimension(12,2) :: eq
       integer, dimension(12) :: ids
       real(real64), dimension(18) :: q
       real(real64), dimension(18) :: expected_f
       real(real64), dimension(18,18) :: expected_m
       real(real64), dimension(18,18) :: expected_k
       real(real64), dimension(12,12) :: local
       real(real64), dimension(12) :: local_f
       real(real64), allocatable, dimension(:) :: f
       real(real64), allocatable, dimension(:) :: threaded_f
       real(real64), allocatable, dimension(:,:) :: dense
       real(real64), allocatable, dimension(:,:) :: threaded_dense
       integer :: i
       integer :: p
       integer :: c

       call mesh%assemble_mass(sparse)
       ok = size(sparse%row_indices) == 1 .and. size(sparse%values) == 0
       call mesh%assemble_stiffness(q(:0), threaded_sparse, nthreads=2)
       ok = ok .and. all(threaded_sparse%row_indices == [1]) .and. size(threaded_sparse%values) == 0
       call mesh%assemble_internal_force(q(:0), threaded_f, nthreads=2)
       ok = ok .and. size(threaded_f) == 0
       mesh%nodes = [fe_node(41, [0.0d0, 0.0d0, 0.0d0]), &
           fe_node(80, [0.0d0, 5.0d0, 0.0d0]), fe_node(25, [2.0d0, 0.0d0, 0.0d0]), &
           fe_node(64, [3.0d0, 0.0d0, 0.0d0])]
       material = linear_elastic_isotropic_material(100.0d0, 0.25d0, 2.0d0)
       beam%material => material
       beam%area = 1.0d0
       beam%j = 0.1d0
       beam%iyy = 0.2d0
       beam%izz = 0.3d0
       allocate(mesh%elements(2))
       beam%node_indices = [1, 3]
       allocate(mesh%elements(1)%element, source=beam)
       beam%node_indices = [4, 3]
       allocate(mesh%elements(2)%element, source=beam)
       call mesh%number_dofs(offsets)
       ok = ok .and. all(offsets == [0, 6, 6, 12, 18])

       eq(:,1) = [(i, i=1,12)]
       eq(:,2) = [(i, i=13,18), (i, i=7,12)]
       q = [(0.001d0*i, i=1,18)]
       expected_m = 0.0d0
       expected_k = 0.0d0
       expected_f = 0.0d0
       do i = 1, 2
           ids = eq(:,i)
           associate (element => mesh%elements(i)%element)
               call element%mass_matrix(local, mesh%nodes)
               expected_m(ids,ids) = expected_m(ids,ids) + local
               call element%stiffness_matrix(q(ids), local, mesh%nodes)
               expected_k(ids,ids) = expected_k(ids,ids) + local
               call element%internal_force_vector(q(ids), local_f, mesh%nodes)
               expected_f(ids) = expected_f(ids) + local_f
           end associate
       end do
       call mesh%assemble_mass(dense)
       ok = ok .and. maxval(abs(dense-expected_m)) < 1.0d-12
       call mesh%assemble_mass(sparse)
       ok = ok .and. maxval(abs(csr_to_dense(sparse)-expected_m)) < 1.0d-12
       call mesh%assemble_stiffness(q, dense)
       ok = ok .and. maxval(abs(dense-expected_k)) < 1.0d-12
       call mesh%assemble_stiffness(q, sparse)
       ok = ok .and. maxval(abs(csr_to_dense(sparse)-expected_k)) < 1.0d-12
       do i = 1, 18
           do p = sparse%row_indices(i)+1, sparse%row_indices(i+1)-1
               c = sparse%column_indices(p)
               ok = ok .and. c > sparse%column_indices(p-1)
           end do
       end do
       call mesh%assemble_internal_force(q, f)
       ok = ok .and. maxval(abs(f-expected_f)) < 1.0d-12

       call mesh%assemble_stiffness(q, threaded_dense, nthreads=2)
       ok = ok .and. all(threaded_dense == dense)
       call mesh%assemble_stiffness(q, threaded_sparse, nthreads=2)
       ok = ok .and. all(threaded_sparse%row_indices == sparse%row_indices)
       if (size(threaded_sparse%values) == size(sparse%values)) then
           ok = ok .and. all(threaded_sparse%column_indices == sparse%column_indices)
           ok = ok .and. all(threaded_sparse%values == sparse%values)
       else
           ok = .false.
       end if
       call mesh%assemble_internal_force(q, threaded_f, nthreads=2)
       ok = ok .and. all(threaded_f == f)
       if (.not.ok) print "(A)", "TEST FAILED: flexible mesh assembly"
   end function

   function test_flexible_mesh_boundary_conditions() result(ok)
       logical :: ok
       type(linear_elastic_isotropic_material), target :: material
       type(euler_bernoulli_beam_element) :: beam
       type(finite_element_mesh) :: mesh
       type(nodal_condition), allocatable, dimension(:) :: conditions
       type(csr_matrix) :: sparse
       integer, dimension(12,2) :: eq
       integer, dimension(12) :: ids
       integer, allocatable, dimension(:) :: free
       real(real64), dimension(18) :: q
       real(real64), dimension(18) :: utarget
       real(real64), dimension(18) :: u
       real(real64), dimension(12,12) :: local
       real(real64), allocatable, dimension(:) :: f
       real(real64), allocatable, dimension(:) :: fr
       real(real64), allocatable, dimension(:) :: fs
       real(real64), allocatable, dimension(:) :: uf
       real(real64), allocatable, dimension(:,:) :: k
       real(real64), dimension(18,18) :: expected_k
       real(real64), dimension(18) :: expected_f
       integer :: i

       ! Same mesh as test_flexible_mesh: node 2 is unused, so there are 18 DOFs
       mesh%nodes = [fe_node(41, [0.0d0, 0.0d0, 0.0d0]), &
           fe_node(80, [0.0d0, 5.0d0, 0.0d0]), fe_node(25, [2.0d0, 0.0d0, 0.0d0]), &
           fe_node(64, [3.0d0, 0.0d0, 0.0d0])]
       material = linear_elastic_isotropic_material(100.0d0, 0.25d0, 2.0d0)
       beam%material => material
       beam%area = 1.0d0
       beam%j = 0.1d0
       beam%iyy = 0.2d0
       beam%izz = 0.3d0
       allocate(mesh%elements(2))
       beam%node_indices = [1, 3]
       allocate(mesh%elements(1)%element, source=beam)
       beam%node_indices = [4, 3]
       allocate(mesh%elements(2)%element, source=beam)

       eq(:,1) = [(i, i=1,12)]
       eq(:,2) = [(i, i=13,18), (i, i=7,12)]
       q = [(0.001d0*i, i=1,18)]
       expected_k = 0.0d0
       expected_f = 0.0d0
       do i = 1, 2
           ids = eq(:,i)
           associate (element => mesh%elements(i)%element)
               call element%stiffness_matrix(q(ids), local, mesh%nodes)
               expected_k(ids,ids) = expected_k(ids,ids) + local
           end associate
       end do

       ! Force conditions: node 3 component 2 is global DOF 8, node 4 component 1 is DOF 13
       conditions = [nodal_condition(3, 2, 2.0d0), nodal_condition(3, 2, 0.5d0), &
           nodal_condition(4, 1, -1.0d0)]
       call apply_force_boundary_conditions(mesh, conditions, f)
       ok = size(f) == 18
       ok = ok .and. abs(f(8) - 2.5d0) < 1.0d-12 .and. abs(f(13) + 1.0d0) < 1.0d-12
       ok = ok .and. abs(sum(abs(f)) - 3.5d0) < 1.0d-12
       call apply_force_boundary_conditions(mesh, conditions(3:3), f)
       ok = ok .and. abs(f(13) + 2.0d0) < 1.0d-12

       ! Displacement conditions on global DOFs 1, 2, and 13
       utarget = [(0.01d0*i, i=1,18)]
       conditions = [nodal_condition(1, 1, utarget(1)), nodal_condition(1, 2, utarget(2)), &
           nodal_condition(4, 1, utarget(13))]
       free = [(i, i=3,12), (i, i=14,18)]
       expected_f = matmul(expected_k, utarget)
       k = expected_k
       fr = expected_f
       call apply_displacement_boundary_conditions(mesh, conditions, k, fr)
       ok = ok .and. size(k,1) == 15 .and. size(fr) == 15
       ok = ok .and. maxval(abs(k - expected_k(free,free))) < 1.0d-12
       ok = ok .and. maxval(abs(fr - matmul(expected_k(free,free), utarget(free)))) < 1.0d-12

       ! Restoring the reduced solution reproduces the prescribed state
       uf = utarget(free)
       u = restore_displacement_boundary_conditions(mesh, conditions, uf)
       ok = ok .and. maxval(abs(u - utarget)) < 1.0d-12

       ! The CSR overload must match the dense reduction
       call mesh%assemble_stiffness(q, sparse)
       fs = expected_f
       call apply_displacement_boundary_conditions(mesh, conditions, sparse, fs)
       ok = ok .and. maxval(abs(csr_to_dense(sparse) - expected_k(free,free))) < 1.0d-12
       ok = ok .and. maxval(abs(fs - matmul(expected_k(free,free), utarget(free)))) < 1.0d-12

       ok = test_labeled_dof_lookup() .and. ok
       if (.not.ok) print "(A)", "TEST FAILED: flexible mesh boundary conditions"
   end function

    function test_labeled_dof_lookup() result(ok)
        logical :: ok
        type(finite_element_mesh) :: mesh
        type(labeled_test_beam) :: beam
        type(degree_of_freedom), allocatable, dimension(:) :: dofs
        type(nodal_condition), dimension(3) :: conditions
        integer, allocatable, dimension(:) :: equations
        real(real64), allocatable, dimension(:) :: f
        integer :: ndof

        mesh%nodes = [fe_node(41, [0.0d0, 0.0d0, 0.0d0]), fe_node(25, [1.0d0, 0.0d0, 0.0d0])]
        beam%node_indices = [1, 2]
        allocate(mesh%elements(1))
        allocate(mesh%elements(1)%element, source=beam)

        ! UY is the first slot and UX is the second in this element.
        conditions = [nodal_condition(2, DOF_UX, 2.0d0), nodal_condition(1, DOF_UY, 1.0d0), &
            nodal_condition(2, DOF_UX, 3.0d0)]
        call mesh%resolve_dofs(conditions%node, conditions%component, dofs, ndof)
        ok = ndof == 12 .and. all(dofs%equation_id == [8, 1, 8])
        ok = ok .and. all(dofs%node_id == [25, 41, 25]) .and. all(dofs%component == conditions%component)
        call resolve_boundary_conditions(mesh, conditions, equations)
        ok = ok .and. all(equations == [8, 1, 8])
        call apply_force_boundary_conditions(mesh, conditions, f)
        ok = ok .and. abs(f(8) - 5.0d0) < 1.0d-12 .and. abs(f(1) - 1.0d0) < 1.0d-12
    end function
end module
