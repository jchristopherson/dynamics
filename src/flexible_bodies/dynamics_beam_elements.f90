module dynamics_flexible_beam_elements
    !! Energy-based three-dimensional flexible beam elements.
    !!
    !! Use a consistent system of length, mass, and time units. Internal
    !! forces are the positive gradient of stored energy and the tangent
    !! stiffness is its Hessian with respect to the element DOFs.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use autodiff
    use dynamics_helper, only : dot_product, cross_product
    use dynamics_finite_element
    use dynamics_flexible_materials
    use dynamics_flexible_body_types, only : fe_node, DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ
    implicit none
    private
    public :: flexible_beam_element
    public :: euler_bernoulli_beam_element
    public :: timoshenko_beam_element

    type, abstract, extends(finite_element) :: flexible_beam_element
        !! Base for beams whose force and tangent derive from strain energy.
        !! Each instance owns a reusable AD workspace; concurrent calls on
        !! the same instance require external synchronization.
        class(flexible_material), pointer :: material => null()
            !! Associated material. The caller must keep the target alive
            !! for as long as the beam uses it.
        real(real64) :: area = 0.0d0
            !! Undeformed cross-sectional area, in units of length squared.
        type(hyperdual_workspace), private :: workspace
            !! Reusable hyper-dual input storage for differentiation.
    contains
        procedure(beam_strain_energy), deferred :: strain_energy
            !! Evaluate the stored energy on hyper-dual element DOFs.
        procedure :: internal_force_vector => beam_internal_force
            !! Evaluate the positive energy gradient.
        procedure :: stiffness_matrix => beam_stiffness
            !! Evaluate the energy Hessian.
    end type

    abstract interface
        function beam_strain_energy(this, u, nodes) result(energy)
            !! Contract for an element's stored strain energy. Implementations
            !! must propagate all hyper-dual components of the input DOFs.
            import flexible_beam_element, hyperdual, fe_node
            class(flexible_beam_element), intent(in) :: this
                !! Beam supplying its geometry and constitutive response.
            type(hyperdual), intent(in), dimension(:) :: u
                !! Element DOFs in the ordering defined by the concrete beam.
            type(fe_node), intent(in), dimension(:) :: nodes
                !! Caller-owned reference node table.
            type(hyperdual) :: energy
                !! Stored energy and its seeded first and mixed derivatives.
        end function
    end interface

    type, abstract, extends(flexible_beam_element) :: two_node_beam_element
        !! Shared layout and section data for spatial beams referenced to
        !! centroidal principal section axes.
        !! The twelve DOFs are ordered [ux, uy, uz, rx, ry, rz] per node.
        !! Translations and rotation vectors are in global axes; rotations
        !! are expressed in radians. Coordinates and section dimensions
        !! describe the undeformed configuration.
        integer, dimension(2) :: node_indices = 0
            !! One-based positions in the supplied node table, in element
            !! DOF order. These are not `fe_node%id` values.
        real(real64), dimension(3) :: reference_y = [0.0d0, 1.0d0, 0.0d0]
            !! Global direction defining local section y after projection
            !! normal to the reference chord; it must not be parallel to it.
        real(real64) :: j = 0.0d0
            !! Positive Saint-Venant torsion constant (length to the fourth).
        real(real64) :: iyy = 0.0d0
            !! Positive second moment of area about local y (length to the fourth).
        real(real64) :: izz = 0.0d0
            !! Positive second moment of area about local z (length to the fourth).
        integer, dimension(12) :: equation_indices = 0
            !! Caller-supplied indices for manual assembly. The mesh numbers
            !! equations from node connectivity instead.
    contains
        procedure :: dof_count => beam_dof_count
        procedure :: dof_indices => beam_dof_indices
        procedure :: node_count => beam_node_count
        procedure :: node_connectivity => beam_node_connectivity
        procedure :: dof_per_node => beam_dof_per_node
        procedure :: dof_components => beam_dof_components
        procedure :: dimensionality => beam_dimensionality
    end type

    type, extends(two_node_beam_element) :: euler_bernoulli_beam_element
        !! Shear-free spatial beam with four-point section energy integration
        !! and constant reference mass.
    contains
        procedure :: mass_matrix => eb_mass
        procedure :: strain_energy => eb_energy
    end type

    type, extends(two_node_beam_element) :: timoshenko_beam_element
        !! Finite-strain spatial beam with linear centerline interpolation
        !! and independently interpolated section rotations. Axial, bending
        !! and torsional energy use four Gauss points; transverse shear uses
        !! the midpoint to relieve thin-beam shear locking.
        real(real64) :: shear_factor_y = 5.0d0/6.0d0
            !! Positive shear correction for displacement along local y.
        real(real64) :: shear_factor_z = 5.0d0/6.0d0
            !! Positive shear correction for displacement along local z.
    contains
        procedure :: mass_matrix => tm_mass
        procedure :: strain_energy => tm_energy
    end type

    type :: beam_context
        !! Read-only beam reference passed through the AD callback.
        class(flexible_beam_element), pointer :: beam
            !! Beam whose energy is being differentiated.
        type(fe_node), pointer, dimension(:) :: nodes
            !! Node table borrowed only for the duration of the Hessian call.
    end type

    interface euler_bernoulli_beam_element
        module procedure :: eb_init
    end interface

    interface timoshenko_beam_element
        module procedure :: tm_init
    end interface

! ------------------------------------------------------------------------------
    real(real64), parameter :: frame_tol = 1.0d-12
    real(real64), parameter, dimension(4) :: points = [ &
        0.06943184420297371d0, &
        0.33000947820757187d0, &
        0.66999052179242813d0, &
        0.93056815579702629d0 &
    ]
    real(real64), parameter, dimension(4) :: weights = [ &
        0.17392742256872693d0, &
        0.32607257743127307d0, &
        0.32607257743127307d0, &
        0.17392742256872693d0 &
    ]

contains
! ------------------------------------------------------------------------------
    pure function beam_dof_components(this) result(components)
        !! Return the global translation and rotation components in local order.
        class(two_node_beam_element), intent(in) :: this
        integer, allocatable, dimension(:) :: components

        components = [DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ]
    end function

! ------------------------------------------------------------------------------
    function eb_init(mat, A, J, Iyy, Izz, n1, n2, ydir) result(rst)
        !! Constructs a new [[euler_bernoulli_beam_element]] instance.
        class(flexible_material), intent(in), pointer :: mat
            !! A pointer to the material to associate with the element.  The
            !! caller must keep the target alive for at least as long as the
            !! lifetime of this object.
        real(real64), intent(in) :: A
            !! The beam cross-sectional area (length squared).
        real(real64), intent(in) :: J
            !! Positive Saint-Venant torsion constant (length to the fourth).
        real(real64), intent(in) :: Iyy
            !! Positive second moment of area about local y (length to the fourth).
        real(real64), intent(in) :: Izz
            !! Positive second moment of area about local z (length to the fourth).
        integer, intent(in) :: n1
            !! Node 1 index in the node table. This value does not necessarily
            !! equal the 'fe_node%id' value.
        integer, intent(in) :: n2
            !! Node 2 index in the node table. This value does not necessarily
            !! equal the 'fe_node%id' value.
        real(real64), intent(in), dimension(3) :: ydir
            !! Global direction defining local section y after projection
            !! normal to the reference chord; it must not be parallel to it.
        type(euler_bernoulli_beam_element) :: rst
            !! The new [[euler_bernoulli_beam_element]] object.

        call init_beam_section(rst, mat, A, J, Iyy, Izz, n1, n2, ydir)
    end function

    function tm_init(mat, A, J, Iyy, Izz, n1, n2, ydir, shear_factor_y, shear_factor_z) result(rst)
        !! Construct a [[timoshenko_beam_element]]. Geometry and section
        !! properties are validated when the element is evaluated.
        class(flexible_material), intent(in), pointer :: mat
            !! Material pointer; its target must outlive the element.
        real(real64), intent(in) :: A
            !! Positive reference cross-sectional area.
        real(real64), intent(in) :: J
            !! Positive Saint-Venant torsion constant.
        real(real64), intent(in) :: Iyy
            !! Positive second moment of area about local y.
        real(real64), intent(in) :: Izz
            !! Positive second moment of area about local z.
        integer, intent(in) :: n1
            !! One-based position of the first node in the node table.
        integer, intent(in) :: n2
            !! One-based position of the second node in the node table.
        real(real64), intent(in), dimension(3) :: ydir
            !! Global direction for local y, projected normal to the chord.
        real(real64), intent(in), optional :: shear_factor_y
            !! Positive, dimensionless local-y shear correction; defaults to 5/6.
        real(real64), intent(in), optional :: shear_factor_z
            !! Positive, dimensionless local-z shear correction; defaults to 5/6.
        type(timoshenko_beam_element) :: rst
            !! Configured element with unset manual equation indices.

        call init_beam_section(rst, mat, A, J, Iyy, Izz, n1, n2, ydir)
        if (present(shear_factor_y)) rst%shear_factor_y = shear_factor_y
        if (present(shear_factor_z)) rst%shear_factor_z = shear_factor_z
    end function

    subroutine init_beam_section(beam, mat, A, J, Iyy, Izz, n1, n2, ydir)
        !! Set the geometry, orientation and material fields shared by the
        !! two beam constructors; borrow the material without copying it.
        class(two_node_beam_element), intent(inout) :: beam
        class(flexible_material), intent(in), pointer :: mat
        real(real64), intent(in) :: A
        real(real64), intent(in) :: J
        real(real64), intent(in) :: Iyy
        real(real64), intent(in) :: Izz
        integer, intent(in) :: n1
        integer, intent(in) :: n2
        real(real64), intent(in), dimension(3) :: ydir

        beam%material => mat
        beam%area = A
        beam%j = J
        beam%iyy = Iyy
        beam%izz = Izz
        beam%node_indices = [n1, n2]
        beam%reference_y = ydir
    end subroutine

! ------------------------------------------------------------------------------
    subroutine beam_internal_force(this, q, f, nodes)
        !! Evaluate the positive strain-energy gradient by hyper-dual AD.
        !! One energy evaluation is required per DOF. Both arrays must have
        !! `dof_count()` entries; mismatches terminate with an error.
        class(flexible_beam_element), intent(inout), target :: this
            !! Beam; its differentiation workspace may be resized.
        real(real64), intent(in), dimension(:) :: q
            !! Element DOFs in the concrete element's ordering.
        real(real64), intent(out), dimension(:) :: f
            !! Positive energy-gradient vector, conjugate to `q`.
        type(fe_node), intent(in), dimension(:) :: nodes
            !! Caller-owned reference node table.

        type(hyperdual) :: energy
        integer :: i, n

        n = this%dof_count()
        if (size(q) /= n .or. size(f) /= n) error stop "Invalid beam force dimensions"
        if (n == 0) return
        if (.not.this%workspace%validate_size(n, 1)) call this%workspace%resize(n, 1)

        do i = 1, n
            this%workspace%xd_pool(i,1) = hyperdual(q(i), 0.0d0, 0.0d0, 0.0d0)
        end do
        do i = 1, n
            this%workspace%xd_pool(i,1)%dx = 1.0d0
            energy = this%strain_energy(this%workspace%xd_pool(1:n,1), nodes)
            f(i) = energy%dx
            this%workspace%xd_pool(i,1)%dx = 0.0d0
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine beam_stiffness(this, q, k, nodes)
        !! Evaluate the symmetric strain-energy Hessian by hyper-dual AD.
        !! Invalid DOF or matrix dimensions terminate with an error.
        class(flexible_beam_element), intent(inout), target :: this
            !! Beam; its differentiation workspace may be resized.
        real(real64), intent(in), dimension(:) :: q
            !! Element DOFs in the concrete element's ordering.
        real(real64), intent(out), dimension(:,:) :: k
            !! `dof_count()`-by-`dof_count()` tangent stiffness matrix.
        type(fe_node), intent(in), target, dimension(:) :: nodes
            !! Caller-owned reference node table.

        type(beam_context) :: ctx
        integer :: n

        n = this%dof_count()
        if (size(q) /= n .or. size(k,1) /= n .or. size(k,2) /= n) &
            error stop "Invalid beam stiffness dimensions"

        ctx%beam => this
        ctx%nodes => nodes
        call hessian(beam_energy_wrapper, q, k, work=this%workspace, args=ctx)
    end subroutine

! ------------------------------------------------------------------------------
    function beam_energy_wrapper(q, args) result(energy)
        !! Adapt the type-bound energy method to the AD scalar callback.
        !! Missing or incompatible context terminates with an error.
        type(hyperdual), intent(in), dimension(:) :: q
            !! Seeded element DOFs.
        class(*), intent(in), optional :: args
            !! A [[beam_context]] containing an associated beam.
        type(hyperdual) :: energy
            !! Beam energy and seeded derivatives.

        if (.not.present(args)) error stop "Missing beam context"

        select type (args)
        type is (beam_context)
            if (.not.associated(args%beam)) error stop "Missing beam"
            energy = args%beam%strain_energy(q, args%nodes)
        class default
            error stop "Invalid beam context"
        end select
    end function

! ------------------------------------------------------------------------------
    pure function beam_dof_count(this) result(n)
        !! Return the total number of element DOFs.
        class(two_node_beam_element), intent(in) :: this
        integer :: n
            !! Twelve DOFs.
        n = 12
    end function

! ------------------------------------------------------------------------------
    pure function beam_dof_indices(this) result(indices)
        !! Copy the caller-supplied equation indices in element DOF order.
        class(two_node_beam_element), intent(in) :: this
            !! Beam containing the index mapping.
        integer, allocatable, dimension(:) :: indices
            !! Allocated array of twelve equation indices.
        indices = this%equation_indices
    end function

! ------------------------------------------------------------------------------
    pure function beam_node_connectivity(this) result(indices)
        !! Return the two node positions in local DOF block order.
        class(two_node_beam_element), intent(in) :: this
        integer, allocatable, dimension(:) :: indices
        indices = this%node_indices
    end function

! ------------------------------------------------------------------------------
    pure function beam_node_count(this) result(n)
        !! Return the number of element nodes.
        class(two_node_beam_element), intent(in) :: this
        integer :: n
            !! Two nodes.
        n = 2
    end function

! ------------------------------------------------------------------------------
    pure function beam_dof_per_node(this) result(n)
        !! Return the number of DOFs at each node.
        class(two_node_beam_element), intent(in) :: this
        integer :: n
            !! Three translations and three rotation-vector components.
        n = 6
    end function

! ------------------------------------------------------------------------------
    pure function beam_dimensionality(this) result(n)
        !! Return the spatial dimension of the element.
        class(two_node_beam_element), intent(in) :: this
        integer :: n
            !! Three spatial dimensions.
        n = 3
    end function

! ------------------------------------------------------------------------------
    subroutine reference_frame(this, nodes, length, frame, chord)
        !! Validate the material association, section and reference geometry,
        !! then construct an orthonormal, right-handed reference frame.
        !! Invalid geometry or nonpositive section properties terminate with
        !! an error.
        class(two_node_beam_element), intent(in) :: this
            !! Beam supplying node indices, orientation, and section properties.
        type(fe_node), intent(in), dimension(:) :: nodes
            !! Caller-owned reference node table.
        real(real64), intent(out) :: length
            !! Positive reference chord length.
        real(real64), intent(out), dimension(3,3) :: frame
            !! Global reference x, y and z unit directors, one per column.
        real(real64), intent(out), optional, dimension(3) :: chord
            !! Unnormalized reference chord from node one to node two.

        real(real64), dimension(3) :: y
        real(real64) :: yn
        integer :: first
        integer :: second

        if (.not.associated(this%material)) error stop "Beam material is not associated"
        if (.not.all(ieee_is_finite([this%area, this%j, this%iyy, this%izz])) .or. &
            min(this%area, this%j, this%iyy, this%izz) <= 0.0d0) &
            error stop "Beam area, J, Iyy, and Izz must be finite and positive"
        first = this%node_indices(1)
        second = this%node_indices(2)
        if (any(this%node_indices < 1) .or. any(this%node_indices > size(nodes))) &
            error stop "Beam node indices are outside the node table"
        if (first == second) error stop "Beam requires two distinct node indices"
        if (.not.all(ieee_is_finite(nodes(first)%coordinates)) .or. &
            .not.all(ieee_is_finite(nodes(second)%coordinates)) .or. &
            .not.all(ieee_is_finite(this%reference_y))) error stop "Invalid beam reference geometry"

        frame(:,1) = nodes(second)%coordinates - nodes(first)%coordinates
        if (present(chord)) chord = frame(:,1)
        length = norm2(frame(:,1))
        if (.not.(length > 0.0d0)) error stop "Beam reference length must be positive"

        frame(:,1) = frame(:,1) / length
        yn = norm2(this%reference_y)
        if (.not.(yn > 0.0d0)) error stop "Beam reference y direction must be nonzero"

        y = this%reference_y / yn
        y = y - dot_product(y, frame(:,1))*frame(:,1)
        yn = norm2(y)
        if (.not.(yn > frame_tol)) error stop "Beam reference y direction is parallel to its axis"

        frame(:,2) = y / yn
        frame(:,3) = cross_product(frame(:,1), frame(:,2))
    end subroutine

! ------------------------------------------------------------------------------
    subroutine eb_mass(this, m, nodes)
        !! Compute the constant 12-by-12 consistent reference mass matrix.
        !! Four-point Gauss integration combines axial and torsional linear
        !! interpolation with Hermite bending interpolation and section
        !! rotary inertia. The polar mass moment uses Iy + Iz, not the
        !! Saint-Venant torsion constant. The material density must be finite
        !! and nonnegative; invalid geometry or matrix shape is an error.
        class(euler_bernoulli_beam_element), intent(inout) :: this
            !! Beam providing its reference geometry and material density.
        real(real64), intent(out), dimension(:,:) :: m
            !! Reference mass matrix in global element DOF order.
        type(fe_node), intent(in), dimension(:) :: nodes
            !! Caller-owned reference node table.

        real(real64), dimension(3,3) :: frame
        real(real64), dimension(6,12) :: b, weighted
        real(real64), dimension(6) :: inertia
        real(real64), dimension(4) :: h, dh
        real(real64) :: length, s
        integer :: i, j, k

        if (size(m,1) /= 12 .or. size(m,2) /= 12) error stop "Invalid beam mass dimensions"
        call reference_frame(this, nodes, length, frame)
        if (.not.ieee_is_finite(this%material%density) .or. this%material%density < 0.0d0) &
            error stop "Invalid beam density"

        inertia = [this%area, this%area, this%area, this%iyy+this%izz, this%iyy, this%izz]
        m = 0.0d0
        do i = 1, 4
            s = points(i)
            b = 0.0d0
            do j = 1, 6
                b(j,j) = 1.0d0 - s
                b(j,j+6) = s
            end do
            h = [1.0d0-3*s**2+2*s**3, length*s*(1-s)**2, 3*s**2-2*s**3, length*s**2*(s-1)]
            dh = [6*s*(s-1)/length, 1-4*s+3*s**2, 6*s*(1-s)/length, 3*s**2-2*s]
            b(2,[2,6,8,12]) = h
            b(3,[3,5,9,11]) = h * [1,-1,1,-1]
            b(5,[3,5,9,11]) = dh * [-1,1,-1,1]
            b(6,[2,6,8,12]) = dh
            do k = 1, 12, 3
                b(:,k:k+2) = matmul(b(:,k:k+2), transpose(frame))
            end do
            weighted = b * spread(inertia, 2, 12)
            m = m + weights(i) * matmul(transpose(b), weighted)
        end do
        m = m * length * this%material%density
    end subroutine

! ------------------------------------------------------------------------------
    function eb_energy(this, u, nodes) result(w)
        !! Integrate the shear-free beam strain energy in the reference
        !! configuration. A cubic Hermite centerline is constructed from the
        !! chord and rotated end tangents; a projected interpolated director
        !! defines the section frame. Axial Green-Lagrange strain is evaluated
        !! on the centroidal axis; the material supplies bending and torsion
        !! energies from the section rotation rates and properties.
        !!
        !! The axial, bending, and torsional material responses are uncoupled.
        !! Zero current chord, singular tangent or section frame, and invalid
        !! reference geometry or DOF count terminate with an error.
        class(euler_bernoulli_beam_element), intent(in) :: this
            !! Beam with an associated material and valid reference section.
        type(hyperdual), intent(in), dimension(:) :: u
            !! Twelve seeded global translations and rotation-vector DOFs.
        type(fe_node), intent(in), dimension(:) :: nodes
            !! Caller-owned reference node table.
        type(hyperdual) :: w
            !! Total stored energy and its seeded derivatives.

        real(real64) :: length, s
        real(real64), dimension(3,3) :: frame
        real(real64), dimension(3) :: reference_chord
        type(hyperdual), dimension(4,2) :: q
        type(hyperdual), dimension(3,2) :: a, b
        type(hyperdual), dimension(3) :: dr, ddr, tangent, dt, guide, dguide, &
            d2, d3, chord
        type(hyperdual) :: scale, speed, projection, guide_norm, kx, ky, kz, &
            axial_strain, section
        integer :: i, j

        if (size(u) /= 12) error stop "Beam requires twelve degrees of freedom"
        call reference_frame(this, nodes, length, frame, reference_chord)
        do i = 1, 3
            chord(i) = reference_chord(i) + (u(i+6) - u(i))
        end do

        scale = dot_product(chord, chord)
        if (.not.(scale%x > 0.0d0)) error stop "Beam current chord is zero"

        scale = sqrt(scale)
        q(:,1) = rotation_quaternion(u(4:6))
        q(:,2) = rotation_quaternion(u(10:12))
        do j = 1, 2
            a(:,j) = rotate_axis(q(:,j), frame(:,1))
            b(:,j) = rotate_axis(q(:,j), frame(:,2))
        end do
        do i = 1, 3
            dguide(i) = (b(i,2) - b(i,1)) / length
        end do

        w = 0.0d0
        do j = 1, 4
            s = points(j)
            do i = 1, 3
                dr(i) = (6*s*(1-s)*chord(i) + scale*((1-4*s+3*s*s)*a(i,1) + (3*s*s-2*s)*a(i,2))) / length
                ddr(i) = (6*(1-2*s)*chord(i) + scale*((6*s-4)*a(i,1) + (6*s-2)*a(i,2))) / length**2
                guide(i) = (1-s)*b(i,1) + s*b(i,2)
            end do
            speed = dot_product(dr, dr)
            if (.not.(speed%x > frame_tol**2)) error stop "Beam tangent is singular"
            speed = sqrt(speed)
            do i = 1, 3
                tangent(i) = dr(i) / speed
            end do
            projection = dot_product(tangent, ddr)
            do i = 1, 3
                dt(i) = (ddr(i) - projection*tangent(i)) / speed
            end do
            projection = dot_product(tangent, guide)
            do i = 1, 3
                d2(i) = guide(i) - projection*tangent(i)
            end do
            guide_norm = dot_product(d2, d2)
            if (.not.(guide_norm%x > frame_tol**2)) error stop "Beam section frame is singular"
            guide_norm = sqrt(guide_norm)
            do i = 1, 3
                d2(i) = d2(i) / guide_norm
            end do
            d3 = cross_product(tangent, d2)
            kx = (dot_product(d3, dguide) - projection*dot_product(d3, dt)) / guide_norm
            ky = -dot_product(d3, dt)
            kz = dot_product(d2, dt)
            axial_strain = 0.5d0*(speed*speed - 1.0d0)
            section = this%area*this%material%plane_stress_energy(axial_strain) + &
                this%material%bending_energy(ky, kz, this%iyy, this%izz) + &
                this%material%torsional_energy(kx, this%j)
            w = w + weights(j)*length*section
        end do
    end function

! ------------------------------------------------------------------------------
    subroutine tm_mass(this, m, nodes)
        !! Consistent reference mass for independently linear translational
        !! and rotational fields. Rotary inertia uses `Iyy + Izz` about the
        !! beam axis, not `J`. Invalid geometry, shape, or a nonfinite or
        !! negative density causes `error stop`.
        class(timoshenko_beam_element), intent(inout) :: this
            !! Beam supplying the section geometry and reference density.
        real(real64), intent(out), dimension(:,:) :: m
            !! Constant 12-by-12 mass matrix in global element DOF order.
        type(fe_node), intent(in), dimension(:) :: nodes
            !! Caller-owned reference node table.

        real(real64), dimension(3,3) :: frame
        real(real64), dimension(6,6) :: inertia
        real(real64), dimension(3) :: rotary
        real(real64) :: length
        integer :: i

        if (size(m,1) /= 12 .or. size(m,2) /= 12) error stop "Invalid beam mass dimensions"
        call reference_frame(this, nodes, length, frame)
        if (.not.ieee_is_finite(this%material%density) .or. this%material%density < 0.0d0) &
            error stop "Invalid beam density"
        rotary = [this%iyy+this%izz, this%iyy, this%izz]
        inertia = 0.0d0
        do i = 1, 3
            inertia(i,i) = this%area
            inertia(4:6,4:6) = inertia(4:6,4:6) + rotary(i)*spread(frame(:,i),2,3)*spread(frame(:,i),1,3)
        end do
        inertia = (this%material%density*length/6.0d0)*inertia
        m(1:6,1:6) = 2.0d0*inertia
        m(7:12,7:12) = 2.0d0*inertia
        m(1:6,7:12) = inertia
        m(7:12,1:6) = inertia
    end subroutine

    function tm_energy(this, u, nodes) result(w)
        !! Integrate the centroidal Green-Lagrange strain, section curvature,
        !! and engineering shear in the rotating section frame, per reference
        !! length. The custom material's 3D law supplies shear energy through
        !! its pure-transverse-shear slice. Invalid geometry, section data,
        !! or DOF count causes `error stop`, as does a relative nodal rotation
        !! too close to pi where the shortest quaternion branch is ambiguous.
        class(timoshenko_beam_element), intent(in) :: this
            !! Beam with an associated material and valid section properties.
        type(hyperdual), intent(in), dimension(:) :: u
            !! Twelve seeded global translations and rotation-vector DOFs.
        type(fe_node), intent(in), dimension(:) :: nodes
            !! Caller-owned reference node table.
        type(hyperdual) :: w
            !! Total stored energy with propagated derivatives.

        real(real64), dimension(3,3) :: frame
        real(real64), dimension(3) :: reference_chord
        real(real64) :: length
        type(hyperdual), dimension(4,2) :: q
        type(hyperdual), dimension(4) :: section_q
        type(hyperdual), dimension(3) :: dr
        type(hyperdual), dimension(3) :: spin
        type(hyperdual), dimension(3) :: curvature
        type(hyperdual), dimension(6) :: shear
        type(hyperdual) :: axial
        type(hyperdual) :: alignment
        type(hyperdual) :: norm_squared
        integer :: i
        integer :: j

        if (size(u) /= 12) error stop "Beam requires twelve degrees of freedom"
        call reference_frame(this, nodes, length, frame, reference_chord)
        if (.not.all(ieee_is_finite([this%shear_factor_y, this%shear_factor_z])) .or. &
            min(this%shear_factor_y, this%shear_factor_z) <= 0.0d0) error stop "Invalid beam shear correction factors"
        do i = 1, 3
            dr(i) = (reference_chord(i) + (u(i+6) - u(i))) / length
        end do
        axial = dot_product(dr, dr)
        if (.not.(axial%x > 0.0d0)) error stop "Beam current chord is zero"
        axial = 0.5d0*(axial - 1.0d0)

        q(:,1) = rotation_quaternion(u(4:6))
        q(:,2) = rotation_quaternion(u(10:12))
        alignment = dot_product(q(:,1), q(:,2))
        if (.not.(abs(alignment%x) > frame_tol)) error stop "Beam relative rotation is too close to pi"
        if (alignment%x < 0.0d0) then
            do i = 1, 4
                q(i,2) = -q(i,2)
            end do
        end if

        ! Spatial spin is 2*vec(q' * conj(q)). For a normalized linear
        ! quaternion blend, the normalization derivative contributes no vector part.
        spin = cross_product(q(2:4,1), q(2:4,2))
        do i = 1, 3
            spin(i) = 2.0d0*(spin(i) + q(1,1)*q(i+1,2) - q(1,2)*q(i+1,1))/length
        end do
        w = this%area*this%material%plane_stress_energy(axial)
        do j = 1, size(points)
            call blend_quaternion(q, points(j), section_q, norm_squared)
            do i = 1, 3
                curvature(i) = dot_product(rotate_axis(section_q, frame(:,i)), spin)/norm_squared
            end do
            w = w + weights(j)*(this%material%bending_energy(curvature(2), curvature(3), this%iyy, this%izz) + &
                this%material%torsional_energy(curvature(1), this%j))
        end do

        ! Selective integration replaces shear strains by their midpoint
        ! values. The material sees the engineering shear entries [xy, zx].
        call blend_quaternion(q, 0.5d0, section_q, norm_squared)
        do i = 1, 6
            shear(i) = 0.0d0
        end do
        shear(4) = sqrt(this%shear_factor_y)*dot_product(rotate_axis(section_q, frame(:,2)), dr)
        shear(6) = sqrt(this%shear_factor_z)*dot_product(rotate_axis(section_q, frame(:,3)), dr)
        w = length*(w + this%area*this%material%energy_density(shear))
    end function

    subroutine blend_quaternion(nodal, s, q, norm_squared)
        !! Normalize the linear blend of sign-aligned nodal quaternions.
        type(hyperdual), intent(in), dimension(4,2) :: nodal
            !! Scalar-first unit quaternions, one per column, in the same hemisphere.
        real(real64), intent(in) :: s
            !! Fractional reference distance from the first node, in [0, 1].
        type(hyperdual), intent(out), dimension(4) :: q
            !! Interpolated unit quaternion with propagated derivatives.
        type(hyperdual), intent(out) :: norm_squared
            !! Squared norm before normalization, used in the curvature formula.
        type(hyperdual) :: norm
        integer :: i

        do i = 1, 4
            q(i) = (1.0d0-s)*nodal(i,1) + s*nodal(i,2)
        end do
        norm_squared = dot_product(q, q)
        norm = sqrt(norm_squared)
        do i = 1, 4
            q(i) = q(i)/norm
        end do
    end subroutine

! ------------------------------------------------------------------------------
    function rotation_quaternion(r) result(q)
        !! Convert a global rotation vector to a scalar-first unit quaternion.
        !! The small-angle series avoids a singular derivative at zero.
        type(hyperdual), intent(in), dimension(3) :: r
            !! Rotation vector in radians, including derivative seeds.
        type(hyperdual), dimension(4) :: q
            !! Quaternion [scalar, x, y, z] with propagated derivatives.

        type(hyperdual) :: a2, half, scale
        integer :: i
        a2 = dot_product(r, r)
        if (a2%x < 1.0d-8) then
            q(1) = 1.0d0 - a2/8.0d0 + a2*a2/384.0d0
            scale = 0.5d0 - a2/48.0d0 + a2*a2/3840.0d0
        else
            half = 0.5d0*sqrt(a2)
            q(1) = cos(half)
            scale = 0.5d0*sin(half)/half
        end if
        do i = 1, 3
            q(i+1) = scale*r(i)
        end do
    end function

! ------------------------------------------------------------------------------
    function rotate_axis(q, a) result(b)
        !! Rotate a real reference director by a hyper-dual unit quaternion.
        type(hyperdual), intent(in), dimension(4) :: q
            !! Scalar-first unit quaternion.
        real(real64), intent(in), dimension(3) :: a
            !! Reference director in global coordinates.
        type(hyperdual), dimension(3) :: b
            !! Rotated director with propagated derivatives.

        type(hyperdual), dimension(3) :: ah, cross, double_cross
        integer :: i

        do i = 1, 3
            ah(i) = a(i)
        end do
        cross = cross_product(q(2:4), ah)
        double_cross = cross_product(q(2:4), cross)
        do i = 1, 3
            b(i) = a(i) + 2.0d0*(q(1)*cross(i) + double_cross(i))
        end do
    end function

! ------------------------------------------------------------------------------
end module
