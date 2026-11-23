module dynamics_flexible_discrete_elements
    !! Follower connectors and centered nodal inertia for flexible meshes.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dynamics_finite_element, only : finite_element
    use dynamics_flexible_body_types, only : fe_node, DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ
    use dynamics_error_handling, only : DYN_INVALID_INPUT_ERROR, DYN_ARRAY_SIZE_ERROR
    use lapack, only : DSYEV
    implicit none
    private
    public :: flexible_spring_element, flexible_damper_element, flexible_point_mass_element

    type, abstract, extends(finite_element) :: discrete_element_base
        integer, dimension(6) :: equation_indices = 0
            !! Optional equation indices for manual assembly; the mesh ignores them.
    contains
        procedure :: dof_count => discrete_dof_count
        procedure :: dof_indices => discrete_dof_indices
        procedure :: dimensionality => spatial_dimension
        procedure :: mass_matrix => zero_mass
        procedure :: stiffness_matrix => zero_stiffness
        procedure :: internal_force_vector => zero_force
    end type

    type, abstract, extends(discrete_element_base) :: axial_connector
        integer, dimension(2) :: node_indices = 0
            !! Distinct one-based positions in the mesh node table.
    contains
        procedure :: node_count => connector_node_count
        procedure :: node_connectivity => connector_nodes
        procedure :: dof_per_node => connector_dof_per_node
        procedure :: dof_components => connector_components
        procedure :: uses_dof_subset => connector_subset
    end type

    type, extends(axial_connector) :: flexible_spring_element
        !! Axial spring following the current chord, stress-free at its reference length.
        real(real64) :: stiffness = 0.0d0
            !! Nonnegative elastic stiffness (force per length).
    contains
        procedure :: internal_force_vector => spring_force
        procedure :: stiffness_matrix => spring_stiffness
    end type

    type, extends(axial_connector) :: flexible_damper_element
        !! Viscous damper following the current chord; no static stiffness.
        real(real64) :: damping_coefficient = 0.0d0
            !! Nonnegative viscous coefficient (force per extension rate).
    contains
        procedure :: add_damping_force => damper_force
    end type

    type, extends(discrete_element_base) :: flexible_point_mass_element
        !! Centered mass and constant rotary inertia at a single six-DOF node.
        integer :: node_index = 0
            !! One-based mesh node-table position.
        real(real64) :: mass = 0.0d0
            !! Nonnegative translational mass.
        real(real64), dimension(3,3) :: inertia = 0.0d0
            !! Symmetric positive semidefinite inertia in global reference axes.
    contains
        procedure :: node_count => mass_node_count
        procedure :: node_connectivity => mass_node
        procedure :: dof_per_node => mass_dof_per_node
        procedure :: dof_components => mass_components
        procedure :: mass_matrix => point_mass_matrix
    end type

    interface flexible_spring_element
        module procedure :: make_spring
    end interface

    interface flexible_damper_element
        module procedure :: make_damper
    end interface

    interface flexible_point_mass_element
        module procedure :: make_point_mass
    end interface

contains

    pure function make_spring(n1, n2, stiffness) result(element)
        !! Construct a follower spring between two node-table positions.
        integer, intent(in) :: n1
        integer, intent(in) :: n2
        real(real64), intent(in) :: stiffness
        type(flexible_spring_element) :: element
        element%node_indices = [n1, n2]
        element%stiffness = stiffness
    end function

    pure function make_damper(n1, n2, damping_coefficient) result(element)
        !! Construct a follower damper between two node-table positions.
        integer, intent(in) :: n1
        integer, intent(in) :: n2
        real(real64), intent(in) :: damping_coefficient
        type(flexible_damper_element) :: element
        element%node_indices = [n1, n2]
        element%damping_coefficient = damping_coefficient
    end function

    pure function make_point_mass(node, mass, inertia) result(element)
        !! Construct a centered nodal mass with optional constant reference inertia.
        integer, intent(in) :: node
        real(real64), intent(in) :: mass
        real(real64), intent(in), optional, dimension(3,3) :: inertia
        type(flexible_point_mass_element) :: element
        element%node_index = node
        element%mass = mass
        if (present(inertia)) element%inertia = inertia
    end function

    pure integer function discrete_dof_count(this) result(n)
        !! Return the total DOF count in node-block order.
        class(discrete_element_base), intent(in) :: this
        n = this%node_count()*this%dof_per_node()
    end function

    pure function discrete_dof_indices(this) result(indices)
        !! Return the optional manual equation indices.
        class(discrete_element_base), intent(in) :: this
        integer, allocatable, dimension(:) :: indices
        indices = this%equation_indices
    end function

    pure integer function spatial_dimension(this) result(n)
        !! Return the spatial dimension.
        class(discrete_element_base), intent(in) :: this
        n = 3
    end function

    pure integer function connector_node_count(this) result(n)
        !! Return the number of attachment nodes.
        class(axial_connector), intent(in) :: this
        n = 2
    end function

    pure function connector_nodes(this) result(indices)
        !! Return the attachment nodes in local DOF-block order.
        class(axial_connector), intent(in) :: this
        integer, allocatable, dimension(:) :: indices
        indices = this%node_indices
    end function

    pure integer function connector_dof_per_node(this) result(n)
        !! Return the number of translation DOFs per attachment.
        class(axial_connector), intent(in) :: this
        n = 3
    end function

    pure function connector_components(this) result(components)
        !! Return the global translation labels in local order.
        class(axial_connector), intent(in) :: this
        integer, allocatable, dimension(:) :: components
        components = [DOF_UX, DOF_UY, DOF_UZ]
    end function

    pure logical function connector_subset(this) result(subset)
        !! Map translations by label without introducing rotations.
        class(axial_connector), intent(in) :: this
        subset = .true.
    end function

    pure integer function mass_node_count(this) result(n)
        !! Return the number of attachment nodes.
        class(flexible_point_mass_element), intent(in) :: this
        n = 1
    end function

    pure function mass_node(this) result(indices)
        !! Return the single attachment node.
        class(flexible_point_mass_element), intent(in) :: this
        integer, allocatable, dimension(:) :: indices
        indices = [this%node_index]
    end function

    pure integer function mass_dof_per_node(this) result(n)
        !! Return the number of translation and rotation DOFs at the node.
        class(flexible_point_mass_element), intent(in) :: this
        n = 6
    end function

    pure function mass_components(this) result(components)
        !! Return the global translation and rotation labels in local order.
        class(flexible_point_mass_element), intent(in) :: this
        integer, allocatable, dimension(:) :: components
        components = [DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ]
    end function

    subroutine zero_mass(this, m, nodes)
        !! Return the zero mass contribution of a connector.
        class(discrete_element_base), intent(inout) :: this
        real(real64), intent(out), dimension(:,:) :: m
        type(fe_node), intent(in), dimension(:) :: nodes
        if (any(shape(m) /= this%dof_count())) error stop DYN_ARRAY_SIZE_ERROR
        m = 0.0d0
    end subroutine

    subroutine zero_stiffness(this, q, k, nodes)
        !! Return the zero static tangent of a damper or nodal mass.
        class(discrete_element_base), intent(inout), target :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: k
        type(fe_node), intent(in), target, dimension(:) :: nodes
        if (size(q) /= this%dof_count() .or. any(shape(k) /= size(q))) error stop DYN_ARRAY_SIZE_ERROR
        k = 0.0d0
    end subroutine

    subroutine zero_force(this, q, f, nodes)
        !! Return the zero displacement-dependent force of a damper or mass.
        class(discrete_element_base), intent(inout), target :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: f
        type(fe_node), intent(in), dimension(:) :: nodes
        if (size(q) /= this%dof_count() .or. size(f) /= size(q)) error stop DYN_ARRAY_SIZE_ERROR
        f = 0.0d0
    end subroutine

    subroutine connector_geometry(this, q, nodes, axis, length, reference_length)
        !! Resolve the current chord and unit axis; a collapsed chord is invalid.
        !! A positive reference length is required when requested by a spring.
        class(axial_connector), intent(in) :: this
        real(real64), intent(in), dimension(:) :: q
        type(fe_node), intent(in), dimension(:) :: nodes
        real(real64), intent(out), dimension(3) :: axis
        real(real64), intent(out) :: length
        real(real64), intent(out), optional :: reference_length
        real(real64), dimension(3) :: chord
        integer :: n1
        integer :: n2

        if (size(q) /= 6) error stop DYN_ARRAY_SIZE_ERROR
        n1 = this%node_indices(1)
        n2 = this%node_indices(2)
        if (min(n1,n2) < 1 .or. max(n1,n2) > size(nodes) .or. n1 == n2) error stop DYN_INVALID_INPUT_ERROR
        chord = nodes(n2)%coordinates - nodes(n1)%coordinates
        if (.not.all(ieee_is_finite(chord))) error stop DYN_INVALID_INPUT_ERROR
        if (present(reference_length)) then
            reference_length = norm2(chord)
            if (.not.ieee_is_finite(reference_length) .or. reference_length <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
        end if
        chord = chord + (q(4:6) - q(1:3))
        if (.not.all(ieee_is_finite(chord))) error stop DYN_INVALID_INPUT_ERROR
        length = norm2(chord)
        if (.not.ieee_is_finite(length) .or. length <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
        axis = chord/length
    end subroutine

    subroutine spring_force(this, q, f, nodes)
        !! Return the gradient of k*(length - reference_length)**2/2.
        class(flexible_spring_element), intent(inout), target :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: f
        type(fe_node), intent(in), dimension(:) :: nodes
        real(real64), dimension(3) :: axis
        real(real64) :: length
        real(real64) :: reference_length

        if (size(f) /= 6) error stop DYN_ARRAY_SIZE_ERROR
        if (.not.ieee_is_finite(this%stiffness) .or. this%stiffness < 0.0d0) error stop DYN_INVALID_INPUT_ERROR
        call connector_geometry(this, q, nodes, axis, length, reference_length)
        f(4:6) = (this%stiffness*(length-reference_length))*axis
        f(1:3) = -f(4:6)
    end subroutine

    subroutine spring_stiffness(this, q, k, nodes)
        !! Differentiate both the axial force magnitude and its rotating axis.
        class(flexible_spring_element), intent(inout), target :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: k
        type(fe_node), intent(in), target, dimension(:) :: nodes
        real(real64), dimension(3,3) :: block
        real(real64), dimension(3) :: axis
        real(real64) :: length
        real(real64) :: reference_length
        real(real64) :: geometric
        integer :: i

        if (any(shape(k) /= 6)) error stop DYN_ARRAY_SIZE_ERROR
        if (.not.ieee_is_finite(this%stiffness) .or. this%stiffness < 0.0d0) error stop DYN_INVALID_INPUT_ERROR
        call connector_geometry(this, q, nodes, axis, length, reference_length)
        geometric = this%stiffness*((length-reference_length)/length)
        block = (this%stiffness-geometric)*spread(axis,2,3)*spread(axis,1,3)
        do i = 1, 3
            block(i,i) = block(i,i) + geometric
        end do
        k(1:3,1:3) = block
        k(4:6,4:6) = block
        k(1:3,4:6) = -block
        k(4:6,1:3) = -block
    end subroutine

    subroutine damper_force(this, q, v, f, nodes)
        !! Add a resisting force pair proportional to the current chord length rate.
        class(flexible_damper_element), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(in), dimension(:) :: v
        real(real64), intent(inout), dimension(:) :: f
        type(fe_node), intent(in), dimension(:) :: nodes
        real(real64), dimension(3) :: axis
        real(real64), dimension(3) :: force
        real(real64) :: length

        if (size(v) /= 6 .or. size(f) /= 6) error stop DYN_ARRAY_SIZE_ERROR
        if (.not.ieee_is_finite(this%damping_coefficient) .or. this%damping_coefficient < 0.0d0) &
            error stop DYN_INVALID_INPUT_ERROR
        call connector_geometry(this, q, nodes, axis, length)
        force = (this%damping_coefficient*dot_product(axis,v(4:6)-v(1:3)))*axis
        f(1:3) = f(1:3) - force
        f(4:6) = f(4:6) + force
    end subroutine

    subroutine point_mass_matrix(this, m, nodes)
        !! Return diag(m*I, J_ref) with no translation-rotation coupling.
        !! Reject nonfinite, nonsymmetric or indefinite inertia and negative mass.
        class(flexible_point_mass_element), intent(inout) :: this
        real(real64), intent(out), dimension(:,:) :: m
        type(fe_node), intent(in), dimension(:) :: nodes
        real(real64), dimension(3,3) :: scaled
        real(real64), dimension(3,3) :: spectrum
        real(real64), dimension(3) :: eigenvalues
        real(real64), dimension(8) :: work
        real(real64) :: scale
        real(real64), parameter :: tol = 16.0d0*epsilon(1.0d0)
        integer :: i
        integer :: info

        if (any(shape(m) /= 6)) error stop DYN_ARRAY_SIZE_ERROR
        if (this%node_index < 1 .or. this%node_index > size(nodes)) error stop DYN_INVALID_INPUT_ERROR
        if (.not.ieee_is_finite(this%mass) .or. this%mass < 0.0d0 .or. &
            .not.all(ieee_is_finite(this%inertia))) error stop DYN_INVALID_INPUT_ERROR
        m = 0.0d0
        do i = 1, 3
            m(i,i) = this%mass
        end do
        scale = maxval(abs(this%inertia))
        if (scale == 0.0d0) return
        scaled = this%inertia/scale
        if (maxval(abs(scaled-transpose(scaled))) > tol) error stop DYN_INVALID_INPUT_ERROR
        scaled = 0.5d0*(scaled+transpose(scaled))
        spectrum = scaled
        call DSYEV('N', 'U', 3, spectrum, 3, eigenvalues, work, size(work), info)
        if (info /= 0) error stop DYN_INVALID_INPUT_ERROR
        if (minval(eigenvalues) < -tol) error stop DYN_INVALID_INPUT_ERROR
        m(4:6,4:6) = scale*scaled
    end subroutine
end module
