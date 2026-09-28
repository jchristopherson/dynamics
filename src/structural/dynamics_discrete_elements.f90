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
module dynamics_discrete_elements
    !! Discrete (lumped-parameter) spring, viscous damper, and point mass
    !! elements for 2D and 3D structural analyses.
    !!
    !! Spring and damper elements connect two nodes and act along a single
    !! unit axis \(\boldsymbol{a}\).  By default, the axis runs from node 1
    !! to node 2; however, a user-defined axis may be supplied, which permits
    !! zero-length (coincident node) elements.  Defining
    !! $$ \boldsymbol{b}=\begin{bmatrix}-\boldsymbol{a}\\
    !! \boldsymbol{a}\end{bmatrix}, $$
    !! the global spring stiffness and damper damping matrices are
    !! $$ K = k\,\boldsymbol{b}\boldsymbol{b}^T, \qquad
    !! C = c\,\boldsymbol{b}\boldsymbol{b}^T. $$
    !!
    !! The point mass element is attached to a single node and contributes
    !! a translational mass matrix
    !! $$ M = m\,I. $$
    !!
    !! Each element carries only translational degrees of freedom (x and y
    !! in 2D; x, y, and z in 3D).  When assembled against nodes carrying
    !! additional degrees of freedom (e.g. beam nodes), the element degrees of
    !! freedom map onto the leading translational degrees of freedom of each
    !! node.
    use iso_fortran_env, only : int32, real64
    use linalg, only : csr_matrix, dense_to_csr
    use dynamics_error_handling
    use dynamics_structural, only : node, element
    implicit none
    private
    public :: discrete_element
    public :: two_node_discrete_element
    public :: spring_element
    public :: spring_element_2d
    public :: spring_element_3d
    public :: damper_element
    public :: damper_element_2d
    public :: damper_element_3d
    public :: mass_element
    public :: mass_element_2d
    public :: mass_element_3d
    public :: assemble_damping_matrix
    public :: assemble_discrete_system

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type, extends(element), abstract :: discrete_element
        !! Defines a discrete (lumped-parameter) element.  Discrete elements
        !! have no spatial extent; therefore, the stiffness, mass, and damping
        !! matrices default to zero and are overridden only where the element
        !! contributes to the system.  The inherited material is unused.
    contains
        procedure, public :: get_node_natural_coordinates => &
            de_get_node_natural_coordinates
        procedure, public :: jacobian => de_jacobian
        procedure, public :: stiffness_matrix => de_stiffness_matrix
        procedure, public :: mass_matrix => de_mass_matrix
        procedure, public :: damping_matrix => de_damping_matrix
        procedure, public :: external_force_vector => de_ext_force_vector
    end type

! ------------------------------------------------------------------------------
    type, extends(discrete_element), abstract :: two_node_discrete_element
        !! Defines a discrete element that connects two nodes and acts along
        !! a single axis.
        type(node) :: node_1
            !! The first node of the element (s = -1).
        type(node) :: node_2
            !! The second node of the element (s = 1).
        real(real64), dimension(3) :: direction = 0.0d0
            !! The unit vector defining the user-supplied element axis.  Only
            !! the leading components corresponding to the element
            !! dimensionality are used.  This value is used only if
            !! use_direction is true.
        logical :: use_direction = .false.
            !! True if the element axis is defined by direction; else, false
            !! if the axis is defined by the vector from node 1 to node 2.
    contains
        procedure, public :: get_node_count => tnde_get_node_count
        procedure, public :: get_dof_per_node => tnde_dof_per_node
        procedure, public :: get_node => tnde_get_node
        procedure, public :: evaluate_shape_function => tnde_shape_function
        procedure, public :: shape_function_matrix => &
            tnde_shape_function_matrix
        procedure, public :: strain_displacement_matrix => &
            tnde_strain_disp_matrix
        procedure, public :: axis => tnde_axis
    end type

! ------------------------------------------------------------------------------
    type, extends(two_node_discrete_element), abstract :: spring_element
        !! Defines a linear, axial spring element.  The element "strain" is
        !! the spring elongation and the element "stress" is the spring force
        !! (positive in tension).
        real(real64) :: stiffness
            !! The spring stiffness.
    contains
        procedure, public :: constitutive_matrix => se_constitutive_matrix
        procedure, public :: stiffness_matrix => se_stiffness_matrix
    end type

! ------------------------------------------------------------------------------
    type, extends(spring_element) :: spring_element_2d
        !! Defines a 2D linear, axial spring element with x and y translations
        !! at each node.
    contains
        procedure, public :: get_dimensionality => se2d_dimensionality
    end type

    interface spring_element_2d
        module procedure :: se2d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(spring_element) :: spring_element_3d
        !! Defines a 3D linear, axial spring element with x, y, and z
        !! translations at each node.
    contains
        procedure, public :: get_dimensionality => se3d_dimensionality
    end type

    interface spring_element_3d
        module procedure :: se3d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(two_node_discrete_element), abstract :: damper_element
        !! Defines a linear, axial viscous damper element.  When evaluated
        !! with nodal velocities, the element "strain" is the elongation rate
        !! and the element "stress" is the damper force (positive in
        !! tension).
        real(real64) :: damping_coefficient
            !! The viscous damping coefficient.
    contains
        procedure, public :: constitutive_matrix => dp_constitutive_matrix
        procedure, public :: damping_matrix => dp_damping_matrix
    end type

! ------------------------------------------------------------------------------
    type, extends(damper_element) :: damper_element_2d
        !! Defines a 2D linear, axial viscous damper element with x and y
        !! translations at each node.
    contains
        procedure, public :: get_dimensionality => dp2d_dimensionality
    end type

    interface damper_element_2d
        module procedure :: dp2d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(damper_element) :: damper_element_3d
        !! Defines a 3D linear, axial viscous damper element with x, y, and z
        !! translations at each node.
    contains
        procedure, public :: get_dimensionality => dp3d_dimensionality
    end type

    interface damper_element_3d
        module procedure :: dp3d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(discrete_element), abstract :: mass_element
        !! Defines a translational point mass element attached to a single
        !! node.
        real(real64) :: mass
            !! The mass.
        type(node) :: node_1
            !! The node to which the mass is attached.
    contains
        procedure, public :: get_node_count => me_get_node_count
        procedure, public :: get_dof_per_node => me_dof_per_node
        procedure, public :: get_node => me_get_node
        procedure, public :: evaluate_shape_function => me_shape_function
        procedure, public :: shape_function_matrix => &
            me_shape_function_matrix
        procedure, public :: strain_displacement_matrix => &
            me_strain_disp_matrix
        procedure, public :: constitutive_matrix => me_constitutive_matrix
        procedure, public :: mass_matrix => me_mass_matrix
    end type

! ------------------------------------------------------------------------------
    type, extends(mass_element) :: mass_element_2d
        !! Defines a 2D translational point mass element with x and y
        !! translations.
    contains
        procedure, public :: get_dimensionality => me2d_dimensionality
    end type

    interface mass_element_2d
        module procedure :: me2d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(mass_element) :: mass_element_3d
        !! Defines a 3D translational point mass element with x, y, and z
        !! translations.
    contains
        procedure, public :: get_dimensionality => me3d_dimensionality
    end type

    interface mass_element_3d
        module procedure :: me3d_init
    end interface

! ******************************************************************************
! OVERLOADED ROUTINES
! ------------------------------------------------------------------------------
    interface assemble_damping_matrix
        module procedure :: assemble_damping_matrix_dense
        module procedure :: assemble_damping_matrix_csr
    end interface

    interface assemble_discrete_system
        module procedure :: assemble_discrete_system_dense
        module procedure :: assemble_discrete_system_csr
    end interface

contains
! ******************************************************************************
! ASSEMBLY ROUTINES
! ------------------------------------------------------------------------------
pure function find_node_global_dof(n, nodes) result(rst)
    !! Finds the first global degree of freedom associated with a node.
    class(node), intent(in) :: n
        !! The node for which to search.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    integer(int32) :: rst
        !! The index of the first global degree of freedom of the node, or
        !! zero if the node is not found in the list.

    ! Local Variables
    integer(int32) :: i, offset

    ! Process
    rst = 0
    offset = 0
    do i = 1, size(nodes)
        if (n%index == nodes(i)%index) then
            rst = offset + 1
            return
        end if
        offset = offset + nodes(i)%dof
    end do
end function

! ------------------------------------------------------------------------------
pure function element_dof_map(gdof, elem, nodes) result(rst)
    !! Maps each local degree of freedom of an element onto its global degree
    !! of freedom.  Local degrees of freedom map onto the leading degrees of
    !! freedom of their respective global nodes.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(discrete_element), intent(in) :: elem
        !! The discrete element.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    integer(int32), allocatable, dimension(:) :: rst
        !! The N-element array of global degree of freedom indices, where N is
        !! the total number of element degrees of freedom.

    ! Local Variables
    integer(int32) :: i, k, npn, nnodes, first

    ! Process
    npn = elem%get_dof_per_node()
    nnodes = elem%get_node_count()
    allocate(rst(npn * nnodes))
    do i = 1, nnodes
        first = find_node_global_dof(elem%get_node(i), nodes)
        if (first == 0) error stop DYN_INVALID_INPUT_ERROR
        do k = 1, npn
            rst((i - 1) * npn + k) = first + k - 1
        end do
    end do
    if (any(rst > gdof)) error stop DYN_INDEX_OUT_OF_RANGE
end function

! ------------------------------------------------------------------------------
pure subroutine add_discrete_contributions(gdof, elements, nodes, m, c, k)
    !! Adds the mass, damping, and stiffness contributions of a collection of
    !! discrete elements to the global system matrices.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(discrete_element), intent(in), dimension(:) :: elements
        !! The discrete elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    real(real64), intent(inout), dimension(:,:) :: m
        !! The gdof-by-gdof global mass matrix to update.
    real(real64), intent(inout), dimension(:,:) :: c
        !! The gdof-by-gdof global damping matrix to update.
    real(real64), intent(inout), dimension(:,:) :: k
        !! The gdof-by-gdof global stiffness matrix to update.

    ! Local Variables
    integer(int32) :: eidx
    integer(int32), allocatable, dimension(:) :: map

    ! Process
    do eidx = 1, size(elements)
        map = element_dof_map(gdof, elements(eidx), nodes)
        m(map, map) = m(map, map) + elements(eidx)%mass_matrix()
        c(map, map) = c(map, map) + elements(eidx)%damping_matrix()
        k(map, map) = k(map, map) + elements(eidx)%stiffness_matrix()
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_damping_matrix_dense(gdof, elements, nodes, c)
    !! Assembles a dense global damping matrix from a collection of discrete
    !! elements.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(discrete_element), intent(in), dimension(:) :: elements
        !! The discrete elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    real(real64), allocatable, intent(out), dimension(:,:) :: c
        !! The assembled gdof-by-gdof global damping matrix.

    ! Local Variables
    integer(int32) :: eidx
    integer(int32), allocatable, dimension(:) :: map

    ! Initialization
    allocate(c(gdof, gdof), source = 0.0d0)

    ! Process
    do eidx = 1, size(elements)
        map = element_dof_map(gdof, elements(eidx), nodes)
        c(map, map) = c(map, map) + elements(eidx)%damping_matrix()
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_damping_matrix_csr(gdof, elements, nodes, c)
    !! Assembles a global damping matrix, in CSR format, from a collection of
    !! discrete elements.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(discrete_element), intent(in), dimension(:) :: elements
        !! The discrete elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    type(csr_matrix), intent(out) :: c
        !! The assembled gdof-by-gdof global damping matrix in CSR format.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: cdense

    ! Process
    call assemble_damping_matrix_dense(gdof, elements, nodes, cdense)
    c = dense_to_csr(cdense)
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_discrete_system_dense(gdof, masses, dampers, springs, &
    nodes, m, c, k)
    !! Assembles the dense global mass, damping, and stiffness matrices of a
    !! system composed of discrete elements.  The mass, damping, and
    !! stiffness contributions of every element in each collection are
    !! assembled; therefore, any discrete element type may be supplied in any
    !! collection.  A zero-sized array may be supplied for any collection
    !! that is not required.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(discrete_element), intent(in), dimension(:) :: masses
        !! The mass elements to assemble.
    class(discrete_element), intent(in), dimension(:) :: dampers
        !! The damper elements to assemble.
    class(discrete_element), intent(in), dimension(:) :: springs
        !! The spring elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    real(real64), allocatable, intent(out), dimension(:,:) :: m
        !! The assembled gdof-by-gdof global mass matrix.
    real(real64), allocatable, intent(out), dimension(:,:) :: c
        !! The assembled gdof-by-gdof global damping matrix.
    real(real64), allocatable, intent(out), dimension(:,:) :: k
        !! The assembled gdof-by-gdof global stiffness matrix.

    ! Initialization
    allocate(m(gdof, gdof), c(gdof, gdof), k(gdof, gdof), source = 0.0d0)

    ! Process
    call add_discrete_contributions(gdof, masses, nodes, m, c, k)
    call add_discrete_contributions(gdof, dampers, nodes, m, c, k)
    call add_discrete_contributions(gdof, springs, nodes, m, c, k)
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_discrete_system_csr(gdof, masses, dampers, springs, &
    nodes, m, c, k)
    !! Assembles the global mass, damping, and stiffness matrices, in CSR
    !! format, of a system composed of discrete elements.  The mass, damping,
    !! and stiffness contributions of every element in each collection are
    !! assembled; therefore, any discrete element type may be supplied in any
    !! collection.  A zero-sized array may be supplied for any collection
    !! that is not required.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(discrete_element), intent(in), dimension(:) :: masses
        !! The mass elements to assemble.
    class(discrete_element), intent(in), dimension(:) :: dampers
        !! The damper elements to assemble.
    class(discrete_element), intent(in), dimension(:) :: springs
        !! The spring elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    type(csr_matrix), intent(out) :: m
        !! The assembled gdof-by-gdof global mass matrix in CSR format.
    type(csr_matrix), intent(out) :: c
        !! The assembled gdof-by-gdof global damping matrix in CSR format.
    type(csr_matrix), intent(out) :: k
        !! The assembled gdof-by-gdof global stiffness matrix in CSR format.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: mdense
    real(real64), allocatable, dimension(:,:) :: cdense
    real(real64), allocatable, dimension(:,:) :: kdense

    ! Process
    call assemble_discrete_system_dense(gdof, masses, dampers, springs, &
        nodes, mdense, cdense, kdense)
    m = dense_to_csr(mdense)
    c = dense_to_csr(cdense)
    k = dense_to_csr(kdense)
end subroutine

! ******************************************************************************
! DISCRETE_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function de_get_node_natural_coordinates(this, i) result(rst)
    !! Returns the natural coordinate of the requested node.  A single-node
    !! element has its node at s = 0; a two-node element has its nodes at
    !! s = -1 and s = 1, respectively.
    class(discrete_element), intent(in) :: this
        !! The discrete_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node.
    real(real64), allocatable, dimension(:) :: rst
        !! The natural coordinate of the node.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = this%get_node_count()
    if (i < 1 .or. i > n) error stop DYN_INDEX_OUT_OF_RANGE
    if (n == 1) then
        allocate(rst(1), source = 0.0d0)
    else
        allocate(rst(1), source = -1.0d0 + 2.0d0 * (i - 1) / (n - 1))
    end if
end function

! ------------------------------------------------------------------------------
pure function de_jacobian(this, s) result(rst)
    !! Returns the element Jacobian.  Discrete elements have no spatial
    !! extent; therefore, a unit 1-by-1 Jacobian is returned.
    class(discrete_element), intent(in) :: this
        !! The discrete_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates at which to evaluate the Jacobian.  This
        !! argument is unused and is present for interface compatibility.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 1-by-1 Jacobian matrix.
    allocate(rst(1,1), source = 1.0d0)
end function

! ------------------------------------------------------------------------------
pure function de_stiffness_matrix(this, rule) result(rst)
    !! Returns the element stiffness matrix.  The default implementation
    !! returns a zero-valued matrix.
    class(discrete_element), intent(in) :: this
        !! The discrete_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused and is present for
        !! interface compatibility.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-N stiffness matrix, where N is the total number of element
        !! degrees of freedom.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = this%get_node_count() * this%get_dof_per_node()
    allocate(rst(n, n), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function de_mass_matrix(this, rule) result(rst)
    !! Returns the element mass matrix.  The default implementation returns a
    !! zero-valued matrix.
    class(discrete_element), intent(in) :: this
        !! The discrete_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused and is present for
        !! interface compatibility.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-N mass matrix, where N is the total number of element
        !! degrees of freedom.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = this%get_node_count() * this%get_dof_per_node()
    allocate(rst(n, n), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function de_damping_matrix(this) result(rst)
    !! Returns the element damping matrix.  The default implementation
    !! returns a zero-valued matrix.
    class(discrete_element), intent(in) :: this
        !! The discrete_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The N-by-N damping matrix, where N is the total number of element
        !! degrees of freedom.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = this%get_node_count() * this%get_dof_per_node()
    allocate(rst(n, n), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function de_ext_force_vector(this, q, rule) result(rst)
    !! Returns the element external force vector.  Discrete elements do not
    !! support distributed loads; therefore, a zero-valued vector is
    !! returned.  Concentrated loads should be applied directly to the global
    !! force vector.
    class(discrete_element), intent(in) :: this
        !! The discrete_element object.
    real(real64), intent(in), dimension(:) :: q
        !! The distributed load vector.  This argument is unused and is
        !! present for interface compatibility.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused and is present for
        !! interface compatibility.
    real(real64), allocatable, dimension(:) :: rst
        !! The N-element force vector, where N is the total number of element
        !! degrees of freedom.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = this%get_node_count() * this%get_dof_per_node()
    allocate(rst(n), source = 0.0d0)
end function

! ******************************************************************************
! TWO_NODE_DISCRETE_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure subroutine initialize_two_node(elem, nd1, nd2, direction)
    !! Initializes the nodes and axis of a [[two_node_discrete_element]].
    class(two_node_discrete_element), intent(inout) :: elem
        !! The two_node_discrete_element object to initialize.
    class(node), intent(in) :: nd1
        !! The first node of the element.
    class(node), intent(in) :: nd2
        !! The second node of the element.
    real(real64), intent(in), optional, dimension(:) :: direction
        !! An optional vector defining the element axis.  The vector must
        !! have as many components as the element dimensionality and must be
        !! non-zero; it is normalized internally.  If not supplied, the axis
        !! is defined by the vector from node 1 to node 2, in which case the
        !! nodes must not be coincident.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: mag
    real(real64), dimension(3) :: v

    ! Process
    n = elem%get_dimensionality()
    elem%node_1 = nd1
    elem%node_2 = nd2
    elem%direction = 0.0d0
    elem%use_direction = .false.
    if (present(direction)) then
        if (size(direction) /= n) error stop DYN_ARRAY_SIZE_ERROR
        mag = norm2(direction)
        if (mag == 0.0d0) error stop DYN_INVALID_INPUT_ERROR
        elem%direction(1:n) = direction / mag
        elem%use_direction = .true.
    else
        v = [nd2%x - nd1%x, nd2%y - nd1%y, nd2%z - nd1%z]
        if (norm2(v(1:n)) == 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
pure function tnde_get_node_count(this) result(rst)
    !! Gets the number of nodes for the element.
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    integer(int32) :: rst
        !! The number of nodes.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function tnde_dof_per_node(this) result(rst)
    !! Gets the number of degrees of freedom per node.  This is equal to the
    !! element dimensionality as only translations are considered.
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    integer(int32) :: rst
        !! The number of DOF per node.
    rst = this%get_dimensionality()
end function

! ------------------------------------------------------------------------------
pure function tnde_get_node(this, i) result(rst)
    !! Gets the requested node from the element.
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node to retrieve.  This value must be either
        !! 1 or 2.
    type(node) :: rst
        !! The requested node.
    select case (i)
    case (1)
        rst = this%node_1
    case (2)
        rst = this%node_2
    case default
        error stop DYN_INDEX_OUT_OF_RANGE
    end select
end function

! ------------------------------------------------------------------------------
pure function tnde_shape_function(this, i, s) result(rst)
    !! Evaluates the i-th linear shape function
    !! $$ N_1 = \frac{1}{2}(1 - s), \qquad N_2 = \frac{1}{2}(1 + s). $$
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    integer(int32), intent(in) :: i
        !! The index of the shape function to evaluate.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the shape function.
    real(real64) :: rst
        !! The value of the i-th shape function at s.
    select case (i)
    case (1)
        rst = 0.5d0 * (1.0d0 - s(1))
    case (2)
        rst = 0.5d0 * (1.0d0 + s(1))
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure function tnde_shape_function_matrix(this, s) result(rst)
    !! Computes the translational displacement interpolation matrix
    !! $$ N = \begin{bmatrix} N_1 I & N_2 I \end{bmatrix}. $$
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The D-by-2D shape function matrix, where D is the element
        !! dimensionality.

    ! Local Variables
    integer(int32) :: k, n
    real(real64) :: n1, n2

    ! Process
    n = this%get_dimensionality()
    n1 = this%evaluate_shape_function(1, s)
    n2 = this%evaluate_shape_function(2, s)
    allocate(rst(n, 2 * n), source = 0.0d0)
    do k = 1, n
        rst(k, k) = n1
        rst(k, k + n) = n2
    end do
end function

! ------------------------------------------------------------------------------
pure function tnde_strain_disp_matrix(this, s) result(rst)
    !! Computes the matrix relating global nodal displacements to the axial
    !! elongation of the element
    !! $$ B = \begin{bmatrix} -\boldsymbol{a}^T & \boldsymbol{a}^T
    !! \end{bmatrix}, $$
    !! where \(\boldsymbol{a}\) is the unit element axis.
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the matrix.  This
        !! argument is unused as the elongation is constant along the element.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 1-by-2D matrix, where D is the element dimensionality.

    ! Local Variables
    integer(int32) :: n
    real(real64), allocatable, dimension(:) :: a

    ! Process
    a = this%axis()
    n = size(a)
    allocate(rst(1, 2 * n))
    rst(1, 1:n) = -a
    rst(1, n + 1:2 * n) = a
end function

! ------------------------------------------------------------------------------
pure function tnde_axis(this) result(rst)
    !! Gets the unit vector defining the element axis in the global
    !! coordinate system.
    class(two_node_discrete_element), intent(in) :: this
        !! The two_node_discrete_element object.
    real(real64), allocatable, dimension(:) :: rst
        !! The D-element unit vector, where D is the element dimensionality.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: mag
    real(real64), dimension(3) :: v

    ! Process
    n = this%get_dimensionality()
    if (this%use_direction) then
        v = this%direction
    else
        v = [this%node_2%x - this%node_1%x, this%node_2%y - this%node_1%y, &
            this%node_2%z - this%node_1%z]
    end if
    mag = norm2(v(1:n))
    if (mag == 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    rst = v(1:n) / mag
end function

! ------------------------------------------------------------------------------
pure function axial_projection(elem) result(rst)
    !! Computes the axial projection matrix \(\boldsymbol{b}\boldsymbol{b}^T\)
    !! where \(\boldsymbol{b}^T = B\) is the element elongation matrix.
    class(two_node_discrete_element), intent(in) :: elem
        !! The two_node_discrete_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 2D-by-2D projection matrix, where D is the element
        !! dimensionality.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: b

    ! Process
    b = elem%strain_displacement_matrix([0.0d0])
    rst = matmul(transpose(b), b)
end function

! ******************************************************************************
! SPRING_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function se_constitutive_matrix(this) result(rst)
    !! Returns the matrix relating spring elongation to spring force.
    class(spring_element), intent(in) :: this
        !! The spring_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 1-by-1 matrix containing the spring stiffness.
    allocate(rst(1,1), source = this%stiffness)
end function

! ------------------------------------------------------------------------------
pure function se_stiffness_matrix(this, rule) result(rst)
    !! Computes the global spring stiffness matrix
    !! $$ K = k\,\boldsymbol{b}\boldsymbol{b}^T, \qquad
    !! \boldsymbol{b}=\begin{bmatrix}-\boldsymbol{a}\\
    !! \boldsymbol{a}\end{bmatrix}. $$
    class(spring_element), intent(in) :: this
        !! The spring_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused as the stiffness
        !! matrix is exact.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 2D-by-2D global stiffness matrix, where D is the element
        !! dimensionality.
    rst = this%stiffness * axial_projection(this)
end function

! ------------------------------------------------------------------------------
pure function se2d_init(k, nd1, nd2, direction) result(rst)
    !! Initializes a new [[spring_element_2d]].
    real(real64), intent(in) :: k
        !! The spring stiffness.
    class(node), intent(in) :: nd1
        !! The first node of the element.
    class(node), intent(in) :: nd2
        !! The second node of the element.
    real(real64), intent(in), optional, dimension(:) :: direction
        !! An optional 2-element vector defining the spring axis.  The vector
        !! must be non-zero and is normalized internally.  If not supplied,
        !! the axis is defined by the vector from node 1 to node 2, in which
        !! case the nodes must not be coincident.
    type(spring_element_2d) :: rst
        !! The new [[spring_element_2d]].
    rst%stiffness = k
    call initialize_two_node(rst, nd1, nd2, direction)
end function

! ------------------------------------------------------------------------------
pure function se2d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(spring_element_2d), intent(in) :: this
        !! The spring_element_2d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function se3d_init(k, nd1, nd2, direction) result(rst)
    !! Initializes a new [[spring_element_3d]].
    real(real64), intent(in) :: k
        !! The spring stiffness.
    class(node), intent(in) :: nd1
        !! The first node of the element.
    class(node), intent(in) :: nd2
        !! The second node of the element.
    real(real64), intent(in), optional, dimension(:) :: direction
        !! An optional 3-element vector defining the spring axis.  The vector
        !! must be non-zero and is normalized internally.  If not supplied,
        !! the axis is defined by the vector from node 1 to node 2, in which
        !! case the nodes must not be coincident.
    type(spring_element_3d) :: rst
        !! The new [[spring_element_3d]].
    rst%stiffness = k
    call initialize_two_node(rst, nd1, nd2, direction)
end function

! ------------------------------------------------------------------------------
pure function se3d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(spring_element_3d), intent(in) :: this
        !! The spring_element_3d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 3
end function

! ******************************************************************************
! DAMPER_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function dp_constitutive_matrix(this) result(rst)
    !! Returns the matrix relating damper elongation rate to damper force.
    class(damper_element), intent(in) :: this
        !! The damper_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 1-by-1 matrix containing the damping coefficient.
    allocate(rst(1,1), source = this%damping_coefficient)
end function

! ------------------------------------------------------------------------------
pure function dp_damping_matrix(this) result(rst)
    !! Computes the global damper damping matrix
    !! $$ C = c\,\boldsymbol{b}\boldsymbol{b}^T, \qquad
    !! \boldsymbol{b}=\begin{bmatrix}-\boldsymbol{a}\\
    !! \boldsymbol{a}\end{bmatrix}. $$
    class(damper_element), intent(in) :: this
        !! The damper_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 2D-by-2D global damping matrix, where D is the element
        !! dimensionality.
    rst = this%damping_coefficient * axial_projection(this)
end function

! ------------------------------------------------------------------------------
pure function dp2d_init(c, nd1, nd2, direction) result(rst)
    !! Initializes a new [[damper_element_2d]].
    real(real64), intent(in) :: c
        !! The viscous damping coefficient.
    class(node), intent(in) :: nd1
        !! The first node of the element.
    class(node), intent(in) :: nd2
        !! The second node of the element.
    real(real64), intent(in), optional, dimension(:) :: direction
        !! An optional 2-element vector defining the damper axis.  The vector
        !! must be non-zero and is normalized internally.  If not supplied,
        !! the axis is defined by the vector from node 1 to node 2, in which
        !! case the nodes must not be coincident.
    type(damper_element_2d) :: rst
        !! The new [[damper_element_2d]].
    rst%damping_coefficient = c
    call initialize_two_node(rst, nd1, nd2, direction)
end function

! ------------------------------------------------------------------------------
pure function dp2d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(damper_element_2d), intent(in) :: this
        !! The damper_element_2d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function dp3d_init(c, nd1, nd2, direction) result(rst)
    !! Initializes a new [[damper_element_3d]].
    real(real64), intent(in) :: c
        !! The viscous damping coefficient.
    class(node), intent(in) :: nd1
        !! The first node of the element.
    class(node), intent(in) :: nd2
        !! The second node of the element.
    real(real64), intent(in), optional, dimension(:) :: direction
        !! An optional 3-element vector defining the damper axis.  The vector
        !! must be non-zero and is normalized internally.  If not supplied,
        !! the axis is defined by the vector from node 1 to node 2, in which
        !! case the nodes must not be coincident.
    type(damper_element_3d) :: rst
        !! The new [[damper_element_3d]].
    rst%damping_coefficient = c
    call initialize_two_node(rst, nd1, nd2, direction)
end function

! ------------------------------------------------------------------------------
pure function dp3d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(damper_element_3d), intent(in) :: this
        !! The damper_element_3d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 3
end function

! ******************************************************************************
! MASS_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function me_get_node_count(this) result(rst)
    !! Gets the number of nodes for the element.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    integer(int32) :: rst
        !! The number of nodes.
    rst = 1
end function

! ------------------------------------------------------------------------------
pure function me_dof_per_node(this) result(rst)
    !! Gets the number of degrees of freedom per node.  This is equal to the
    !! element dimensionality as only translations are considered.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    integer(int32) :: rst
        !! The number of DOF per node.
    rst = this%get_dimensionality()
end function

! ------------------------------------------------------------------------------
pure function me_get_node(this, i) result(rst)
    !! Gets the requested node from the element.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node to retrieve.  This value must be 1.
    type(node) :: rst
        !! The requested node.
    if (i /= 1) error stop DYN_INDEX_OUT_OF_RANGE
    rst = this%node_1
end function

! ------------------------------------------------------------------------------
pure function me_shape_function(this, i, s) result(rst)
    !! Evaluates the i-th shape function.  The single shape function of a
    !! point mass is unity.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    integer(int32), intent(in) :: i
        !! The index of the shape function to evaluate.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the shape function.
    real(real64) :: rst
        !! The value of the i-th shape function at s.
    if (i == 1) then
        rst = 1.0d0
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure function me_shape_function_matrix(this, s) result(rst)
    !! Computes the displacement interpolation matrix, which is the D-by-D
    !! identity matrix, where D is the element dimensionality.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The D-by-D shape function matrix.

    ! Local Variables
    integer(int32) :: k, n

    ! Process
    n = this%get_dimensionality()
    allocate(rst(n, n), source = 0.0d0)
    do k = 1, n
        rst(k, k) = this%evaluate_shape_function(1, s)
    end do
end function

! ------------------------------------------------------------------------------
pure function me_strain_disp_matrix(this, s) result(rst)
    !! Computes the strain-displacement matrix.  A point mass does not
    !! deform; therefore, a zero-valued matrix is returned.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the matrix.  This
        !! argument is unused.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 1-by-D zero-valued matrix, where D is the element
        !! dimensionality.
    allocate(rst(1, this%get_dimensionality()), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function me_constitutive_matrix(this) result(rst)
    !! Returns the constitutive matrix.  A point mass does not deform;
    !! therefore, a zero-valued matrix is returned.
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 1-by-1 zero-valued matrix.
    allocate(rst(1,1), source = 0.0d0)
end function

! ------------------------------------------------------------------------------
pure function me_mass_matrix(this, rule) result(rst)
    !! Computes the element mass matrix
    !! $$ M = m\,I. $$
    class(mass_element), intent(in) :: this
        !! The mass_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused as the mass matrix
        !! is exact.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The D-by-D mass matrix, where D is the element dimensionality.

    ! Local Variables
    integer(int32) :: k, n

    ! Process
    n = this%get_dimensionality()
    allocate(rst(n, n), source = 0.0d0)
    do k = 1, n
        rst(k, k) = this%mass
    end do
end function

! ------------------------------------------------------------------------------
pure function me2d_init(m, nd) result(rst)
    !! Initializes a new [[mass_element_2d]].
    real(real64), intent(in) :: m
        !! The mass.
    class(node), intent(in) :: nd
        !! The node to which the mass is attached.
    type(mass_element_2d) :: rst
        !! The new [[mass_element_2d]].
    rst%mass = m
    rst%node_1 = nd
end function

! ------------------------------------------------------------------------------
pure function me2d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(mass_element_2d), intent(in) :: this
        !! The mass_element_2d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function me3d_init(m, nd) result(rst)
    !! Initializes a new [[mass_element_3d]].
    real(real64), intent(in) :: m
        !! The mass.
    class(node), intent(in) :: nd
        !! The node to which the mass is attached.
    type(mass_element_3d) :: rst
        !! The new [[mass_element_3d]].
    rst%mass = m
    rst%node_1 = nd
end function

! ------------------------------------------------------------------------------
pure function me3d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(mass_element_3d), intent(in) :: this
        !! The mass_element_3d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 3
end function

! ------------------------------------------------------------------------------
end module
