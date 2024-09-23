! Shape Functions:
! 2D Line: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp1.html#shp2dlinerdof
! 3D Line: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp2.html#shp3d2node

module dynamics_structural
    use iso_fortran_env
    use linalg, only : csr_matrix, create_csr_matrix, sort
    use ferror
    use dynamics_error_handling
    use dynamics_rotation
    implicit none
    private
    public :: DYN_ONE_POINT_INTEGRATION_RULE
    public :: DYN_TWO_POINT_INTEGRATION_RULE
    public :: DYN_THREE_POINT_INTEGRATION_RULE
    public :: DYN_FOUR_POINT_INTEGRATION_RULE
    public :: node
    public :: material
    public :: element
    public :: line_element
    public :: beam_element_2d
    public :: shape_function_derivative
    public :: shape_function_second_derivative
    public :: create_connectivity_matrix
    public :: apply_boundary_conditions
    public :: apply_displacement_constraint
    public :: restore_constrained_values
    public :: point
    public :: beam_element_3d

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: DYN_ONE_POINT_INTEGRATION_RULE = 1
        !! Defines a single-point integration rule.
    integer(int32), parameter :: DYN_TWO_POINT_INTEGRATION_RULE = 2
        !! Defines a two-point integration rule.
    integer(int32), parameter :: DYN_THREE_POINT_INTEGRATION_RULE = 3
        !! Defines a three-point integration rule.
    integer(int32), parameter :: DYN_FOUR_POINT_INTEGRATION_RULE = 4
        !! Defines a four-point integration rule.

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type :: point
        !! Defines a point in 3D, Cartesian space.
        real(real64) :: x
            !! The x-coordinate.
        real(real64) :: y
            !! The y-coordinate.
        real(real64) :: z
            !! The z-coordinate.
    end type

! ------------------------------------------------------------------------------
    type, extends(point) :: node
        !! Defines a node.
        integer(int32) :: index
            !! The global index of the node.
        integer(int32) :: dof
            !! The number of degrees of freeedom associated with this node.
    end type

! ------------------------------------------------------------------------------
    type :: material
        !! Defines a linear-elastic-isotropic material.
        real(real64) :: density
            !! The density of the material.
        real(real64) :: modulus
            !! The modulus of elasticity of the material.
        real(real64) :: poissons_ratio
            !! The Poisson's ratio of the material.
    end type

! ------------------------------------------------------------------------------
    type, abstract :: element
        !! Defines an element.
        type(material) :: material
            !! The material.
    contains
        procedure(element_query), deferred, public, pass :: get_dimensionality
        procedure(element_query), deferred, public, pass :: get_node_count
        procedure(element_get_node), deferred, public, pass :: get_node
        procedure(element_query), deferred, public, pass :: get_dof_per_node
        procedure(element_shape_function), deferred, public, pass :: &
            evaluate_shape_function
        procedure(element_matrix_function), deferred, public, pass :: &
            shape_function_matrix
        procedure(element_matrix_function), deferred, public, pass :: &
            strain_displacement_matrix
        procedure(element_const_matrix_function), deferred, public, &
            pass :: constitutive_matrix
        procedure(element_matrix_function), deferred, public, pass :: &
            jacobian
        procedure, public :: stiffness_matrix => e_stiffness_matrix
        procedure, public :: mass_matrix => e_mass_matrix
        procedure, public :: external_force_vector => e_ext_force_vector
    end type

! ------------------------------------------------------------------------------
    type, extends(element), abstract :: line_element
        !! Defines a line element type.
        real(real64) :: area
            !! The element cross-sectional area.
    contains
        procedure(line_element_get_terminal), deferred, public, pass :: &
            get_terminal_nodes
        procedure(line_element_const_matrix_function), deferred, public, &
            pass :: rotation_matrix
        procedure, public :: length => le_length
        procedure, public :: stiffness_matrix => le_stiffness_matrix
        procedure, public :: mass_matrix => le_mass_matrix
        procedure, public :: external_force_vector => le_ext_force_vector
    end type

! ------------------------------------------------------------------------------
    type, extends(line_element) :: beam_element_2d
        !! Defines a two-dimensional Bernoulli-Euler beam element.
        real(real64) :: moment_of_inertia
            !! The beam moment of inertia (second moment of area).
        type(node) :: node_1
            !! The first node of the element (s = -1).
        type(node) :: node_2
            !! The second node of the element (s = 1).
    contains
        procedure, public :: get_dimensionality => b2d_dimensionality
        procedure, public :: get_node_count => b2d_get_node_count
        procedure, public :: get_dof_per_node => b2d_dof_per_node
        procedure, public :: get_node => b2d_get_node
        procedure, public :: get_terminal_nodes => b2d_terminal_nodes
        procedure, public :: evaluate_shape_function => b2d_shape_function
        procedure, public :: shape_function_matrix => b2d_shape_function_matrix_2d
        procedure, public :: strain_displacement_matrix => b2d_strain_disp_matrix_2d
        procedure, public :: constitutive_matrix => b2d_constitutive_matrix
        procedure, public :: jacobian => b2d_jacobian
        procedure, public :: rotation_matrix => b2d_rotation_matrix
        procedure, public :: stiffness_matrix => b2d_stiffness_matrix
        procedure, public :: mass_matrix => b2d_mass_matrix
    end type

! ------------------------------------------------------------------------------
    type, extends(line_element) :: beam_element_3d
        !! Defines a three-dimensional Bernoulli-Euler beam element.
        real(real64) :: Ixx
            !! The beam moment of inertia about the element x-axis.
        real(real64) :: Iyy
            !! The beam moment of inertia about the element y-axis.
        real(real64) :: Izz
            !! The beam moment of inertia about the element z-axis.
        type(node) :: node_1
            !! The first node of the element (s = -1).
        type(node) :: node_2
            !! The second node of the element (s = 1).
        type(point) :: orientation_point
            !! A point used to determine the orientation of the beam in 3D
            !! space.  The orientation point is measured relative to the first
            !! node in the element.  Specifically, the element z axis is assumed
            !! to be defined by the location of this point relative to the
            !! location of node 1.
    contains
        procedure, public :: get_dimensionality => b3d_dimensionality
        procedure, public :: get_node_count => b3d_get_node_count
        procedure, public :: get_dof_per_node => b3d_dof_per_node
        procedure, public :: get_node => b3d_get_node
        procedure, public :: get_terminal_nodes => b3d_terminal_nodes
        procedure, public :: evaluate_shape_function => b3d_shape_function
        procedure, public :: shape_function_matrix => b3d_shape_function_matrix_3d
        procedure, public :: strain_displacement_matrix => b3d_strain_disp_matrix_3d
        procedure, public :: constitutive_matrix => b3d_constitutive_matrix
        procedure, public :: jacobian => b3d_jacobian
        procedure, public :: rotation_matrix => b3d_rotation_matrix
        procedure, public :: stiffness_matrix => b3d_stiffness_matrix
        procedure, public :: mass_matrix => b3d_mass_matrix
    end type

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        pure function element_query(this) result(rst)
            !! Defines the signature of a function performing a query on an
            !! integer-valued property of a element type.
            use iso_fortran_env, only : int32
            import element
            class(element), intent(in) :: this
                !! The element object.
            integer(int32) :: rst
                !! The resulting value.
        end function

        pure function element_get_node(this, i) result(rst)
            !! Defines the signature of a function for retrieving the requested
            !! node from the element.
            use iso_fortran_env, only : int32
            import node
            import element
            class(element), intent(in) :: this
                !! The element object.
            integer(int32), intent(in) :: i
                !! The local index of the node to retrieve.
            type(node) :: rst
                !! The node.
        end function

        pure function element_matrix_function(this, s) result(rst)
            !! Defines the signature of a routine for returning a matrix
            !! associated with the element.
            use iso_fortran_env, only : real64
            import element
            class(element), intent(in) :: this
                !! The element object.
            real(real64), intent(in), dimension(:) :: s
                !! The value of the natural coordinates at which the matrix
                !! should be evaluated.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The resulting matrix.
        end function

        pure function element_const_matrix_function(this) result(rst)
            !! Defines the signature of a routine for returning a matrix
            !! associated with the element.
            use iso_fortran_env, only : real64
            import element
            class(element), intent(in) :: this
                !! The element object.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The resulting matrix.
        end function

        pure function element_shape_function(this, i, s) result(rst)
            !! Defines the signature of a routine for computing the value of
            !! the i-th element shape function at natural coordinate.
            use iso_fortran_env, only : int32, real64
            import element
            class(element), intent(in) :: this
                !! The element object.
            integer(int32), intent(in) :: i
                !! The index of the shape function to evaluate.
            real(real64), intent(in), dimension(:) :: s
                !! The value of the natural coordinates at which to evaluate
                !! the shape function.
            real(real64) :: rst
                !! The value of the i-th shape function at s.
        end function
        pure subroutine line_element_get_terminal(this, i1, i2)
            !! Defines the signature of a routine for returning the terminal
            !! node numbers.
            use iso_fortran_env, only : int32
            import line_element
            class(line_element), intent(in) :: this
                !! The line_element object.
            integer(int32), intent(out) :: i1
                !! The index of the node at the head of the element.
            integer(int32), intent(out) :: i2
                !! The index of the node at the tail of the element.
        end subroutine

        pure function line_element_const_matrix_function(this) result(rst)
            !! Defines the signature of a routine for returning a matrix
            !! associated with the line_element.
            use iso_fortran_env, only : real64
            import line_element
            class(line_element), intent(in) :: this
                !! The line_element object.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The resulting matrix.
        end function

        pure function integrand(elem, s) result(rst)
            !! Defines the signature of a function containing an integrand.
            use iso_fortran_env, only : real64
            import element
            class(element), intent(in) :: elem
                !! The element object.
            real(real64), intent(in), dimension(:) :: s
                !! The natural coordinate at which to evaluate the integrand.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The result.
        end function
    end interface

! ******************************************************************************
! OVERLOADED ROUTINES
! ------------------------------------------------------------------------------
    interface apply_boundary_conditions
        module procedure :: apply_boundary_conditions_mtx
        module procedure :: apply_boundary_conditions_vec
    end interface

contains
! ******************************************************************************
! DIFFERENTIATION ROUTINES
! ------------------------------------------------------------------------------
pure function shape_function_derivative(index, elem, s, i) result(rst)
    !! Computes the derivative of the shape function with respect to the natural
    !! coordinate specified.
    integer(int32), intent(in) :: index
        !! The index of the shape function to evaluate.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the derivative.
    integer(int32), intent(in) :: i
        !! The index of the natural coordinate to with which the derivative is
        !! to be computed.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    real(real64) :: na, nb, h(size(s))

    ! Initialization
    h = 0.0d0
    h(i) = sqrt(epsilon(na))

    ! Process
    na = elem%evaluate_shape_function(index, s + h)
    nb = elem%evaluate_shape_function(index, s - h)
    rst = (na - nb) / (2.0d0 * h(i))
end function

! ------------------------------------------------------------------------------
pure function shape_function_second_derivative(index, elem, s, i) result(rst)
    !! Computes the second derivative of the shape function with respect to the
    !! natural coordinate specified.
    integer(int32), intent(in) :: index
        !! The index of the shape function to evaluate.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the derivative.
    integer(int32), intent(in) :: i
        !! The index of the natural coordinate to with which the derivative is
        !! to be computed.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    real(real64) :: na, nb, nc, h(size(s))

    ! Initialization
    h = 0.0d0
    h(i) = (epsilon(na))**0.25d0
    na = elem%evaluate_shape_function(index, s + h)
    nb = elem%evaluate_shape_function(index, s)
    nc = elem%evaluate_shape_function(index, s - h)
    rst = (na - 2.0d0 * nb + nc) / (h(i)**2)
end function

! ******************************************************************************
! INTEGRATION
! ------------------------------------------------------------------------------
pure function get_model_parameters(rule) result(rst)
    !! Gets the requested integration model parameters.
    integer(int32), intent(in) :: rule
        !! The integration rule.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integration parameters.

    ! Local Variables
    real(real64) :: x, w1, w2, pt1, pt2

    ! Process
    select case (rule)
    case (DYN_ONE_POINT_INTEGRATION_RULE)
        allocate(rst(2, 2))
        rst = reshape([0.0d0, 0.0d0, 2.0d0, 2.0d0], [2, 2])
    case (DYN_TWO_POINT_INTEGRATION_RULE)
        allocate(rst(2, 2))
        x = sqrt(3.0d0) / 3.0d0
        rst = reshape([-x, x, 1.0d0, 1.0d0], [2, 2])
    case (DYN_THREE_POINT_INTEGRATION_RULE)
        allocate(rst(3, 2))
        x = sqrt(3.0d0 / 5.0d0)
        rst = reshape([0.0d0, -x, x, w1, w2, w2], [3, 2])
    case default ! Four Point Rule
        allocate(rst(4, 2))
        pt1 = sqrt((3.0d0 / 7.0d0) - (2.0d0 / 7.0d0) * sqrt(6.0d0 / 5.0d0))
        pt1 = sqrt((3.0d0 / 7.0d0) + (2.0d0 / 7.0d0) * sqrt(6.0d0 / 5.0d0))
        w1 = (1.8d1 + sqrt(3.0d1)) / 3.6d1
        w2 = (1.8d1 - sqrt(3.0d1)) / 3.6d1
        rst = reshape([-pt1, pt1, -pt2, pt2, w1, w1, w2, w2], [4, 2])
    end select
end function

! ------------------------------------------------------------------------------
pure function integrate_1d(fcn, elem, rule) result(rst)
    !! Computes the integral of the specified integrand given an element and an
    !! integration rule.
    procedure(integrand) :: fcn
        !! The integrand.
    class(element), intent(in) :: elem
        !! The element object.
    integer(int32), intent(in) :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
    real(real64), allocatable, dimension(:,:) :: rst
        !! The result of the integration.

    ! Local Variables
    integer(int32) :: i
    real(real64), allocatable, dimension(:,:) :: s

    ! Process
    s = get_model_parameters(rule)
    rst = s(1,2) * fcn(elem, [s(1,1)])
    do i = 2, size(s, 1)
        rst = rst + s(i,2) * fcn(elem, [s(i,1)])
    end do
end function

! ------------------------------------------------------------------------------
pure function integrate(fcn, elem, rule) result(rst)
    !! Computes the integral of the specified integrand given an element and an
    !! integration rule.
    procedure(integrand) :: fcn
        !! The integrand.
    class(element), intent(in) :: elem
        !! The element object.
    integer(int32), intent(in) :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
    real(real64), allocatable, dimension(:,:) :: rst
        !! The result of the integration.

    ! Process
    select type (elem)
    class is (line_element)
        rst = integrate_1d(fcn, elem, rule)
    end select
end function

! ******************************************************************************
! ASSEMBLY ROUTINES
! ------------------------------------------------------------------------------
pure function find_global_dof(n, nodes) result(rst)
    !! Finds the index of the global DOF node in a list of nodes.
    class(node), intent(in) :: n
        !! The node for which to search.
    class(node), intent(in), dimension(:) :: nodes
        !! The list of nodes
    integer(int32) :: rst
        !! The requested index.

    ! Local Variables
    integer(int32) :: i

    ! Process
    rst = 0
    do i = 1, size(nodes)
        if (n%index == nodes(i)%index) then
            rst = rst + 1
            exit
        end if
        rst = rst + nodes(i)%dof
    end do
end function

! ------------------------------------------------------------------------------
function create_connectivity_matrix(gdof, e, nodes, err) result(rst)
    !! Creates a connectivity matrix for the element.
    integer(int32), intent(in) :: gdof
        !! The number of global degrees of freedom.
    class(element), intent(in) :: e
        !! The element.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    class(errors), intent(inout), optional, target :: err
        !! An optional error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: i, j, col, nnodes, nnz, row, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nnodes = e%get_node_count()
    nnz = e%get_dof_per_node() * nnodes
    allocate(rst(nnz, gdof), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call report_memory_error("create_connectivity_matrix", flag, errmgr)
        return
    end if

    ! Process
    row = 0
    do j = 1, nnodes
        col = find_global_dof(e%get_node(j), nodes)
        do i = 1, e%get_dof_per_node()
            row = row + 1
            rst(row, col) = 1.0d0
            col = col + 1
        end do
    end do
end function

! ------------------------------------------------------------------------------
function apply_boundary_conditions_mtx(gdof, x, err) result(rst)
    !! Applies boundary conditions to a matrix by removal of the appropriate
    !! rows and columns.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    real(real64), intent(in), dimension(:,:) :: x
        !! The matrix to constrain.
    class(errors), intent(inout), optional, target :: err
        !! An optional error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The altered matrix.

    ! Local Variables
    integer(int32) :: i, j, ii, m, n, nbc, mnew, flag
    integer(int32), allocatable, dimension(:) :: indices
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = size(x, 1)
    n = size(x, 2)
    nbc = size(gdof)
    mnew = m - nbc

    ! Input Checking
    if (m /= n) then
        call report_nonsquare_matrix_error("apply_boundary_conditions_mtx", &
            "x", m, n, errmgr)
        return
    end if
    if (mnew < 1) then
        call report_overconstraint_error("apply_boundary_conditions_mtx", &
            errmgr)
        return
    end if
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > m) then
            call report_array_index_out_of_bounds_error( &
                "apply_boundary_conditions_mtx", "gdof", gdof(i), m, errmgr)
            return
        end if
    end do

    ! Memory Allocation
    allocate(rst(mnew, mnew), stat = flag)
    if (flag == 0) allocate(indices(m - nbc), stat = flag)
    if (flag /= 0) then
        call report_memory_error("apply_boundary_conditions_mtx", flag, errmgr)
        return
    end if

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) then
            call report_nonmonotonic_array_error(&
                "apply_boundary_conditions_mtx", "gdof", i, errmgr)
            return
        end if
    end do

    ! Process
    ii = 1
    j = 0
    do i = 1, m
        if (gdof(ii) /= i) then
            j = j + 1
            indices(j) = i
        else
            ii = ii + 1
            if (ii > nbc) ii = nbc
        end if
    end do

    ! Now, we only need store the rows and columns stored in indices
    rst = x(indices,indices)
end function

! ------------------------------------------------------------------------------
function apply_boundary_conditions_vec(gdof, x, err) result(rst)
    !! Applies boundary conditions to a vector by removal of the appropriate
    !! items.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    real(real64), intent(in), dimension(:) :: x
        !! The vector to constrain.
    class(errors), intent(inout), optional, target :: err
        !! An optional error handling object.
    real(real64), allocatable, dimension(:) :: rst
        !! The altered vector.

    ! Local Variables
    integer(int32) :: i, j, ii, n, nbc, nnew, flag
    integer(int32), allocatable, dimension(:) :: indices
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    nbc = size(gdof)
    nnew = n - nbc

    ! Input Checking
    if (nnew < 1) then
        call report_overconstraint_error("apply_boundary_conditions_vec", &
            errmgr)
        return
    end if
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > n) then
            call report_array_index_out_of_bounds_error( &
                "apply_boundary_conditions_vec", "gdof", gdof(i), n, errmgr)
            return
        end if
    end do

    ! Memory Allocation
    allocate(rst(nnew), stat = flag)
    if (flag == 0) allocate(indices(n - nbc), stat = flag)
    if (flag /= 0) then
        call report_memory_error("apply_boundary_conditions_vec", flag, errmgr)
        return
    end if

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) then
            call report_nonmonotonic_array_error( &
                "apply_boundary_conditions_vec", "gdof", i, errmgr)
            return
        end if
    end do

    ! Process
    ii = 1
    j = 0
    do i = 1, n
        if (gdof(ii) /= i) then
            j = j + 1
            indices(j) = i
        else
            ii = ii + 1
            if (ii > nbc) ii = nbc
        end if
    end do

    ! Now, just store the appropriate items in the output vector
    rst = x(indices)
end function

! ------------------------------------------------------------------------------
function restore_constrained_values(gdof, x, err) result(rst)
    !! Restores the constrained degrees-of-freedom from the boundary conditions
    !! applied by apply_boundary_conditions.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    real(real64), intent(in), dimension(:) :: x
        !! The constrained vector.
    class(errors), intent(inout), optional, target :: err
        !! An optional error handling object.
    real(real64), allocatable, dimension(:) :: rst
        !! The altered vector.

    ! Local Variables
    integer(int32) ::i, j, ii, n, nbc, nnew, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    nbc = size(gdof)
    nnew = n + nbc

    ! Input Checking
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > nnew) then
            call report_array_index_out_of_bounds_error( &
                "restore_constrained_values", "gdof", gdof(i), nnew, errmgr)
            return
        end if
    end do

    ! Memory Allocation
    allocate(rst(nnew), source = 0.0d0, stat = flag)
    if (flag /= 0) then
        call report_memory_error("restore_constrained_values", flag, errmgr)
        return
    end if

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) then
            call report_nonmonotonic_array_error( &
                "restore_constrained_values", "gdof", i, errmgr)
            return
        end if
    end do

    ! Process
    ii = 1
    j = 0
    do i = 1, nnew
        if (i == gdof(ii)) then
            ii = ii + 1
            if (ii > nbc) ii = nbc
        else
            j = j + 1
            rst(i) = x(j)
        end if
    end do
end function

! ------------------------------------------------------------------------------
! REF: https://www.sciencedirect.com/topics/engineering/prescribed-displacement-boundary-condition
subroutine apply_displacement_constraint(dof, val, k, f)
    !! Applies a displacement constraint to the specified degree of freedom.
    integer(int32), intent(in) :: dof
        !! The global degree-of-freedom to which the constraint should be
        !! applied.
    real(real64), intent(in) :: val
        !! The value of the displacement constraint.
    real(real64), intent(inout), dimension(:,:) :: k
        !! The stiffness matrix to which the constraint should be applied.
    real(real64), intent(inout), dimension(:) :: f
        !! The external force vector to which the constraint should be applied.

    ! Wipe out the rows in the matrix and place a value of 1 on the diagonal
    k(dof,:) = 0.0d0
    k(dof,dof) = 1.0d0

    ! Update the external force vector
    f(dof) = val
end subroutine

! ******************************************************************************
! ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function e_stiffness_matrix(this, rule) result(rst)
    !! Computes the stiffness matrix for the element.
    class(element), intent(in) :: this
        !! The element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: r

    ! Initialization
    if (present(rule)) then
        r = rule
    else
        r = DYN_TWO_POINT_INTEGRATION_RULE
    end if

    ! Process
    rst = integrate(element_stiffness_integrand, this, r)
end function

! ----------
pure function element_stiffness_integrand(elem, s) result(rst)
    !! The integrand function for computing the stiffness matrix of an element.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate vector at which to evaluate the integrand.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integrand.

    ! Local Variables
    real(real64) :: jdet
    real(real64), allocatable, dimension(:,:) :: b, bt, d, x, jac

    ! Process
    b = elem%strain_displacement_matrix(s)
    bt = transpose(b)
    d = elem%constitutive_matrix()
    jac = elem%jacobian(s)
    jdet = det(jac)
    x = matmul(d, b)
    rst = jdet * matmul(bt, x)
end function

! ------------------------------------------------------------------------------
pure function e_mass_matrix(this, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(element), intent(in) :: this
        !! The element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: r

    ! Initialization
    if (present(rule)) then
        r = rule
    else
        r = DYN_TWO_POINT_INTEGRATION_RULE
    end if

    ! Process
    rst = integrate(element_mass_integrand, this, r)
end function

! ----------
pure function element_mass_integrand(elem, s) result(rst)
    !! The integrand function for computing the mass matrix of an element.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate vector at which to evaluate the integrand.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integrand.

    ! Local Variables
    real(real64) :: jdet
    real(real64), allocatable, dimension(:,:) :: N, Nt, jac

    ! Process
    N = elem%shape_function_matrix(s)
    Nt = transpose(N)
    jac = elem%jacobian(s)
    jdet = det(jac)
    rst = elem%material%density * jdet * matmul(Nt, N)
end function

! ------------------------------------------------------------------------------
pure function e_ext_force_vector(this, q, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(element), intent(in) :: this
        !! The element object.
    real(real64), intent(in), dimension(:) :: q
        !! The surface traction forces vector or body force vector.  
        !! For instance, a 2D problem this vector would look like [qx, qy]**T.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    integer(int32) :: r

    ! Initialization
    if (present(rule)) then
    r = rule
    else
    r = DYN_TWO_POINT_INTEGRATION_RULE
    end if

    ! Process
    rst = matmul( &
        integrate(element_ext_force_integrand, this, r), &
        q &
    )
end function

! ----------
pure function element_ext_force_integrand(elem, s) result(rst)
    !! The integrand function for computing the external force vector of an 
    !! element.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate vector at which to evaluate the integrand.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integrand.

    ! Local Variables
    real(real64) :: jdet
    real(real64), allocatable, dimension(:,:) :: Nt, jac

    ! Process
    Nt = transpose(elem%shape_function_matrix(s))
    jac = elem%jacobian(s)
    jdet = det(jac)
    rst = jdet * Nt
end function

! ******************************************************************************
! LINE_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function le_length(this) result(rst)
    !! Computes the length of the line_element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    real(real64) :: rst
        !! The length of the line element.

    ! Local Variables
    real(real64) :: dx, dy, dz
    integer(int32) :: i1, i2
    type(node) :: n1, n2

    ! Process
    call this%get_terminal_nodes(i1, i2)
    n1 = this%get_node(i1)
    n2 = this%get_node(i2)
    dx = n2%x - n1%x
    dy = n2%y - n1%y
    dz = n2%z - n1%z
    rst = sqrt(dx**2 + dy**2 + dz**2)
end function


! ------------------------------------------------------------------------------
pure function le_stiffness_matrix(this, rule) result(rst)
    !! Computes the stiffness matrix for the element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Compute the stiffness matrix and apply the rotation transformation
    rst = e_stiffness_matrix(this, rule)
    rst = matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
pure function le_mass_matrix(this, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Compute the mass matrix and apply the rotation transformation
    rst = e_mass_matrix(this, rule)
    rst = this%area * matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
pure function le_ext_force_vector(this, q, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    real(real64), intent(in), dimension(:) :: q
        !! The surface traction forces vector or body force vector.  
        !! For instance, a 2D problem this vector would look like [qx, qy]**T.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: T

    ! Compute the rotation matrix
    T = this%rotation_matrix()

    ! Compute the force vector
    rst = e_ext_force_vector(this, q, rule)
    rst = matmul(T, rst)
end function

! ******************************************************************************
! BEAM_2D ROUTINES
! ------------------------------------------------------------------------------
pure function b2d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function b2d_get_node_count(this) result(rst)
    !! Gets the number of nodes for the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32) :: rst
        !! The number of nodes.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function b2d_dof_per_node(this) result(rst)
    !! Gets the number of degrees of freedom per node.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32) :: rst
        !! The number of DOF per node.
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function b2d_get_node(this, i) result(rst)
    !! Gets the requested node from the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32), intent(in) :: i
        !! The local index of the node to retrieve.
    type(node) :: rst
        !! The requested node.

    if (i == 1) then
        rst = this%node_1
    else
        rst = this%node_2
    end if
end function

! ------------------------------------------------------------------------------
pure subroutine b2d_terminal_nodes(this, i1, i2)
    !! Gets the terminal node numbers for the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32), intent(out) :: i1
        !! The index of the node at the head of the element.
    integer(int32), intent(out) :: i2
        !! The index of the node at the tail of the element.
    i1 = 1
    i2 = 2
end subroutine

! ------------------------------------------------------------------------------
pure function b2d_shape_function(this, i, s) result(rst)
    !! Evaluates the i-th shape function at natural coordinate s.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32), intent(in) :: i
        !! The index of the shape function to evaluate.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate
        !! the shape function.
    real(real64) :: rst
        !! The value of the i-th shape function at s.

    ! Local Variables
    real(real64) :: l

    ! Process
    select case (i)
    case (1)
        rst = 0.5d0 * (1.0d0 - s(1))
    case (2)
        rst = 0.25d0 * (1.0d0 - s(1))**2 * (2.0d0 + s(1))
    case (3)
        rst = 0.25d0 * (1.0d0 - s(1))**2 * (1.0d0 + s(1))
    case (4)
        rst = 0.5d0 * (1.0d0 + s(1))
    case (5)
        rst = 0.25d0 * (1.0d0 + s(1))**2 * (2.0d0 - s(1))
    case (6)
        rst = 0.25d0 * (1.0d0 + s(1))**2 * (s(1) - 1.0d0)
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure function b2d_shape_function_matrix_2d(this, s) result(rst)
    !! Computes the shape function matrix for a beam element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate the shape
        !! functions.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The shape function matrix.

    ! Local Variables
    real(real64) :: n1, n2, n3, n4, n5, n6, l

    ! Initialization
    allocate(rst(2, 6), source = 0.0d0)
    l = this%length()

    ! Process
    n1 = this%evaluate_shape_function(1, s)
    n2 = this%evaluate_shape_function(2, s)
    n3 = 0.5d0 * l * this%evaluate_shape_function(3, s)
    n4 = this%evaluate_shape_function(4, s)
    n5 = this%evaluate_shape_function(5, s)
    n6 = 0.5d0 * l * this%evaluate_shape_function(6, s)
    
    rst(1,1) = n1
    rst(2,2) = n2
    rst(2,3) = n3
    rst(1,4) = n4
    rst(2,5) = n5
    rst(2,6) = n6
end function

! ------------------------------------------------------------------------------
pure function b2d_strain_disp_matrix_2d(this, s) result(rst)
    !! Computes the strain-displacement matrix for a 2D beam element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The strain-displacement matrix.

    ! Local Variables
    real(real64) :: l, dsdx, dn1ds, dn2ds, dn3ds, dn4ds, dn5ds, dn6ds
    
    ! Initialization
    allocate(rst(2, 6), source = 0.0d0)

    ! Process
    l = this%length()
    dsdx = 2.0d0 / l    ! s = 2 * x / L - 1, so ds/dx = 2 / L
    dn1ds = shape_function_derivative(1, this, s, 1)
    dn2ds = shape_function_second_derivative(2, this, s, 1)
    dn3ds = shape_function_second_derivative(3, this, s, 1)
    dn4ds = shape_function_derivative(4, this, s, 1)
    dn5ds = shape_function_second_derivative(5, this, s, 1)
    dn6ds = shape_function_second_derivative(6, this, s, 1)
    rst(1,1) = dn1ds * dsdx
    rst(2,2) = dn2ds * dsdx**2
    rst(2,3) = 0.5d0 * l * dn3ds * dsdx**2
    rst(1,4) = dn4ds * dsdx
    rst(2,5) = dn5ds * dsdx**2
    rst(2,6) = 0.5d0 * l * dn6ds * dsdx**2
end function

! ------------------------------------------------------------------------------
pure function b2d_constitutive_matrix(this) result(rst)
    !! Computes the constitutive matrix for the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Process
    allocate(rst(2,2), source = 0.0d0)
    rst(1,1) = this%area * this%material%modulus
    rst(2,2) = this%moment_of_inertia * this%material%modulus
end function

! ------------------------------------------------------------------------------
pure function b2d_jacobian(this, s) result(rst)
    !! Computes the Jacobian matrix for a 2D beam element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The Jacobian matrix.

    rst = reshape([0.5d0 * this%length()], [1, 1])
end function

! ------------------------------------------------------------------------------
pure function b2d_rotation_matrix(this) result(rst)
    !! Computes the rotation matrix for the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting 6-by-6 rotation matrix.

    ! Local Variables
    real(real64) :: theta, ct, st
    type(node) :: n1, n2

    ! Process
    allocate(rst(6, 6), source = 0.0d0)
    n1 = this%get_node(1)
    n2 = this%get_node(2)
    theta = atan2(n2%y - n1%y, n2%x - n1%x)
    ct = cos(theta)
    st = sin(theta)

    rst(1,1) = ct
    rst(2,1) = st
    rst(1,2) = -st
    rst(2,2) = ct
    rst(3,3) = 1.0d0
    rst(4,4) = ct
    rst(5,4) = st
    rst(4,5) = -st
    rst(5,5) = ct
    rst(6,6) = 1.0d0
end function

! ------------------------------------------------------------------------------
pure function b2d_stiffness_matrix(this, rule) result(rst)
    !! Computes the stiffness matrix for the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64) :: A, E, I, L
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Initialization
    A = this%area
    E = this%material%modulus
    I = this%moment_of_inertia
    L = this%length()

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Construct the stiffness matrix
    allocate(rst(6, 6), source = 0.0d0)
    rst(1,1) = A * E / L
    rst(2,2) = 12.0d0 * E * I / (L**3)
    rst(3,3) = 4.0d0 * E * I / L
    rst(2,3) = 6.0d0 * E * I / (L**2)
    rst(3,2) = rst(2,3)
    rst(1,4) = -rst(1,1)
    rst(4,1) = rst(1,4)
    rst(2,5) = -rst(2,2)
    rst(5,2) = rst(2,5)
    rst(2,6) = rst(2,3)
    rst(6,2) = rst(2,6)
    rst(3,5) = -rst(2,3)
    rst(5,3) = rst(3,5)
    rst(3,6) = 2.0d0 * E * I / L
    rst(6,3) = rst(3,6)
    rst(4:6,4:6) = rst(1:3,1:3)
    rst(5,6) = -rst(2,3)
    rst(6,5) = rst(5,6)

    ! Apply the transformation
    rst = matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
pure function b2d_mass_matrix(this, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(beam_element_2d), intent(in) :: this
        !! The beam_element_2d object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64) :: rho, A, L, f
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Initialization
    rho = this%material%density
    A = this%area
    L = this%length()
    f = rho * A * L / 4.2d2

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Construct the mass matrix
    allocate(rst(6, 6), source = 0.0d0)
    rst(1,1) = 1.4d2 * f
    rst(4,1) = 7.0d1 * f
    rst(2,2) = 1.56d2 * f
    rst(3,2) = 2.2d1 * L * f
    rst(5,2) = 5.4d1 * f
    rst(6,2) = -1.3d1 * L * f
    rst(2,3) = rst(3,2)
    rst(3,3) = 4.0d0 * L**2 * f
    rst(5,3) = 1.3d1 * L * f
    rst(6,3) = -3.0d0 * L**2 * f
    rst(1,4) = rst(4,1)
    rst(4,4) = rst(1,1)
    rst(2,5) = rst(5,2)
    rst(3,5) = rst(5,3)
    rst(5,5) = rst(2,2)
    rst(6,5) = -2.2d1 * L * f
    rst(2,6) = rst(6,2)
    rst(3,6) = rst(6,3)
    rst(5,6) = rst(6,5)
    rst(6,6) = rst(3,3)

    ! Apply the transformation
    rst = matmul(Tt, matmul(rst, T))
end function

! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
pure function det_1(x) result(rst)
    ! Determinant of a 1-by-1 matrix.
    real(real64), intent(in), dimension(:,:) :: x
    real(real64) :: rst
    rst = x(1,1)
end function

! ------------------------------------------------------------------------------
pure function det_2(x) result(rst)
    ! Determinant of a 2-by-2 matrix.
    real(real64), intent(in), dimension(:,:) :: x
    real(real64) :: rst

    rst = x(1,1) * x(2,2) - x(1,2) * x(2,1)
end function

! ------------------------------------------------------------------------------
pure function det_3(x) result(rst)
    ! Determinant of a 3-by-3 matrix.
    real(real64), intent(in), dimension(:,:) :: x
    real(real64) :: rst
    rst = x(1,1) * (x(2,2) * x(3,3) - x(2,3) * x(3,2)) - &
        x(1,2) * (x(2,1) * x(3,3) - x(2,3) * x(3,1)) + &
        x(1,3) * (x(2,1) * x(3,2) - x(2,2) * x(3,1))
end function

! ------------------------------------------------------------------------------
pure function det(x) result(rst)
    !! Computes the determinant of a matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The matrix on which to operate.
    real(real64) :: rst
        !! The determinant.

    select case (size(x, 1))
    case (1)
        rst = det_1(x)
    case (2)
        rst = det_2(x)
    case (3)
        rst = det_3(x)
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
! REF: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
pure function normal_vector_to_line(pt1, pt2, pt) result(rst)
    !! Computes the normal vector to a line defined by pt1 and pt2 assuming
    !! some point (pt) not on the line.
    class(point), intent(in) :: pt1
        !! The origin point of the line segment.
    class(point), intent(in) :: pt2
        !! The termination point of the line segment.
    class(point), intent(in) :: pt
        !! A point, not on the line.
    real(real64) :: rst(3)
        !! The resulting normal vector (unit length).

    ! Local Variables
    real(real64) :: a(3), p(3), n(3), amp(3)

    ! Initialization
    a = [pt1%x, pt1%y, pt1%z]
    n = [pt2%x, pt2%y, pt2%z] - a
    p = [pt%x, pt%y, pt%z]
    amp = a - p
    rst = amp - dot_product(amp, n) * n
    rst = rst / norm2(rst)
end function


! ******************************************************************************
! BEAM_3D ROUTINES
! ------------------------------------------------------------------------------
pure function b3d_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32) :: rst
        !! The dimensionality.
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function b3d_get_node_count(this) result(rst)
    !! Gets the number of nodes for the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32) :: rst
        !! The number of nodes.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function b3d_dof_per_node(this) result(rst)
    !! Gets the number of degrees of freedom per node.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32) :: rst
        !! The number of DOF per node.
    rst = 6
end function

! ------------------------------------------------------------------------------
pure function b3d_get_node(this, i) result(rst)
    !! Gets the requested node from the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32), intent(in) :: i
        !! The local index of the node to retrieve.
    type(node) :: rst
        !! The requested node.

    if (i == 1) then
        rst = this%node_1
    else
        rst = this%node_2
    end if
end function

! ------------------------------------------------------------------------------
pure subroutine b3d_terminal_nodes(this, i1, i2)
    !! Gets the terminal node numbers for the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32), intent(out) :: i1
        !! The index of the node at the head of the element.
    integer(int32), intent(out) :: i2
        !! The index of the node at the tail of the element.
    i1 = 1
    i2 = 2
end subroutine

! ------------------------------------------------------------------------------
pure function b3d_shape_function(this, i, s) result(rst)
    !! Evaluates the i-th shape function at natural coordinate s.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32), intent(in) :: i
        !! The index of the shape function to evaluate.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate
        !! the shape function.
    real(real64) :: rst
        !! The value of the i-th shape function at s.

    ! Local Variables
    real(real64) :: l

    ! Process
    select case (i)
    case (1)
        rst = 0.5d0 * (1.0d0 - s(1))
    case (2)
        rst = 0.25d0 * (1.0d0 - s(1))**2 * (2.0d0 + s(1))
    case (3)
        rst = 0.25d0 * (1.0d0 - s(1))**2 * (2.0d0 + s(1))
    case (4)
        rst = 0.5d0 * (1.0d0 - s(1))
    case (5)
        rst = 0.25d0 * (1.0d0 - s(1)**2) * (1.0d0 - s(1))
    case (6)
        rst = 0.25d0 * (1.0d0 - s(1)**2) * (1.0d0 - s(1))
    case (7)
        rst = 0.5d0 * (1.0d0 + s(1))
    case (8)
        rst = 0.25d0 * (1.0d0 + s(1))**2 * (2.0d0 - s(1))
    case (9)
        rst = 0.25d0 * (1.0d0 + s(1))**2 * (2.0d0 - s(1))
    case (10)
        rst = 0.5d0 * (1.0d0 + s(1))
    case (11)
        rst = 0.25d0 * (1.0d0 - s(1)**2) * (1.0d0 + s(1))
    case (12)
        rst = 0.25d0 * (1.0d0 - s(1)**2) * (1.0d0 + s(1))
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure function b3d_shape_function_matrix_3d(this, s) result(rst)
    !! Computes the shape function matrix for a beam element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate the shape
        !! functions.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The shape function matrix.

    ! Local Variables
    real(real64) :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, l

    ! Initialization
    allocate(rst(4, 12), source = 0.0d0)
    l = this%length()

    ! Process
    n1 = this%evaluate_shape_function(1, s)
    n2 = this%evaluate_shape_function(2, s)
    n3 = this%evaluate_shape_function(3, s)
    n4 = this%evaluate_shape_function(4, s)
    n5 = -0.5d0 * l * this%evaluate_shape_function(5, s)
    n6 = 0.5d0 * l * this%evaluate_shape_function(6, s)
    n7 = this%evaluate_shape_function(7, s)
    n8 = this%evaluate_shape_function(8, s)
    n9 = this%evaluate_shape_function(9, s)
    n10 = this%evaluate_shape_function(10, s)
    n11 = 0.5d0 * l * this%evaluate_shape_function(11, s)
    n12 = -0.5d0 * l * this%evaluate_shape_function(12, s)

    rst(1,1) = n1
    rst(2,2) = n2
    rst(3,3) = n3
    rst(4,4) = n4
    rst(3,5) = n5
    rst(2,6) = n6

    rst(1,7) = n7
    rst(2,8) = n8
    rst(3,9) = n9
    rst(4,10) = n10
    rst(3,11) = n11
    rst(2,12) = n12
end function

! ------------------------------------------------------------------------------
pure function b3d_strain_disp_matrix_3d(this, s) result(rst)
    !! Computes the strain-displacement matrix for a 3D beam element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The strain-displacement matrix.

    ! Local Variables
    real(real64) :: l, dsdx, dn1ds, dn2ds, dn3ds, dn4ds, dn5ds, dn6ds, &
        dn7ds, dn8ds, dn9ds, dn10ds, dn11ds, dn12ds
    
    ! Initialization
    allocate(rst(4, 12), source = 0.0d0)

    ! Process
    l = this%length()
    dsdx = 2.0d0 / l    ! s = 2 * x / L - 1, so ds/dx = 2 / L
    dn1ds = shape_function_derivative(1, this, s, 1)
    dn2ds = shape_function_second_derivative(2, this, s, 1)
    dn3ds = shape_function_second_derivative(3, this, s, 1)
    dn4ds = shape_function_derivative(4, this, s, 1)
    dn5ds = -0.5d0 * l * shape_function_second_derivative(5, this, s, 1)
    dn6ds = 0.5d0 * l * shape_function_second_derivative(6, this, s, 1)
    dn7ds = shape_function_derivative(7, this, s, 1)
    dn8ds = shape_function_second_derivative(8, this, s, 1)
    dn9ds = shape_function_second_derivative(9, this, s, 1)
    dn10ds = shape_function_derivative(10, this, s, 1)
    dn11ds = 0.5d0 * l * shape_function_second_derivative(11, this, s, 1)
    dn12ds = -0.5d0 * l * shape_function_second_derivative(12, this, s, 1)

    ! Build the matrix
    rst(1,1) = dn1ds * dsdx
    rst(2,2) = dn2ds * dsdx**2
    rst(3,3) = dn3ds * dsdx**2
    rst(4,4) = dn4ds * dsdx
    rst(3,5) = dn5ds * dsdx**2
    rst(2,6) = dn6ds * dsdx**2
    rst(1,7) = dn7ds * dsdx
    rst(2,8) = dn8ds * dsdx**2
    rst(3,9) = dn9ds * dsdx**2
    rst(4,10) = dn10ds * dsdx
    rst(3,11) = dn11ds * dsdx**2
    rst(2,12) = dn12ds * dsdx**2
end function
! ------------------------------------------------------------------------------
pure function b3d_constitutive_matrix(this) result(rst)
    !! Computes the constitutive matrix for the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Process
    allocate(rst(4,4), source = 0.0d0)
    rst(1,1) = this%area * this%material%modulus
    rst(2,2) = this%Izz * this%material%modulus
    rst(3,3) = this%Iyy * this%material%modulus
    rst(4,4) = this%Ixx * this%material%modulus / &
        (2.0d0 * (1.0d0 + this%material%poissons_ratio))
end function

! ------------------------------------------------------------------------------
pure function b3d_jacobian(this, s) result(rst)
    !! Computes the Jacobian matrix for a 3D beam element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    real(real64), intent(in), dimension(:) :: s
        !! The value of the natural coordinate at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The Jacobian matrix.

    rst = reshape([0.5d0 * this%length()], [1, 1])
end function

! ------------------------------------------------------------------------------
pure function b3d_rotation_matrix(this) result(rst)
    !! Computes the rotation matrix for the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting 12-by-12 rotation matrix.

    ! Local Variables
    real(real64) :: i(3), j(3), k(3)

    ! Define the unit vectors
    i = [ &
        this%node_2%x - this%node_1%x, &
        this%node_2%y - this%node_1%y, &
        this%node_2%z - this%node_1%z &
    ]
    i = i / norm2(i)
    k = normal_vector_to_line(this%node_1, this%node_2, this%orientation_point)
    j = [ &
        k(2) * i(3) - k(3) * i(2), &
        k(3) * i(1) - k(1) * i(3), &
        k(1) * i(2) - k(2) * i(1) &
    ]

    ! Construct the matrix
    allocate(rst(12, 12), source = 0.0d0)
    rst(1:3,1:3) = rotate(i, j, k)
    rst(4:6,4:6) = rst(1:3,1:3)
    rst(7:9,7:9) = rst(1:3,1:3)
    rst(10:12,10:12) = rst(1:3,1:3)
end function

! ------------------------------------------------------------------------------
! https://www.researchgate.net/publication/352816965_3-D_Beam_Finite_Element_Programming_-A_Practical_Guide_Part_1
! https://homes.civil.aau.dk/jc/FemteSemester/Beams3D.pdf
! https://www.sesamx.io/blog/beam_finite_element/
! https://www.brown.edu/Departments/Engineering/Courses/En2340/Projects/Projects_2015/Wenqiang_Fan.pdf
! https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp2.html#shp3d2node
pure function b3d_stiffness_matrix(this, rule) result(rst)
    !! Computes the stiffness matrix for the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64) :: A, E, Iyy, Izz, Jxx, G, L
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Initialization
    A = this%area
    Jxx = this%Ixx
    Iyy = this%Iyy
    Izz = this%Izz
    E = this%material%modulus
    G = E / (2.0d0 * (1.0d0 + this%material%poissons_ratio))
    L = this%length()

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Construct the stiffness matrix
    allocate(rst(12, 12), source = 0.0d0)
    rst(1,1) = A * E / L
    rst(7,1) = -rst(1,1)
    rst(2,2) = 1.2d1 * E * Izz / L**3
    rst(6,2) = 6.0d0 * E * Izz / L**2
    rst(8,2) = -rst(2,2)
    rst(12,2) = rst(6,2)
    rst(3,3) = 1.2d1 * E * Iyy / L**3
    rst(5,3) = -6.0d0 * E * Iyy / L**2
    rst(9,3) = -rst(3,3)
    rst(11,3) = rst(5,3)
    rst(4,4) = G * Jxx / L
    rst(10,4) = -rst(4,4)
    rst(3,5) = rst(5,3)
    rst(5,5) = 4.0d0 * E * Iyy / L
    rst(9,5) = -rst(11,3)
    rst(11,5) = 2.0d0 * E * Iyy / L
    rst(2,6) = rst(6,2)
    rst(6,6) = 4.0d0 * E * Izz / L
    rst(8,6) = -rst(12,2)
    rst(12,6) = 2.0d0 * E * Izz / L
    rst(1,7) = rst(7,1)
    rst(7,7) = rst(1,1)
    rst(2,8) = rst(8,2)
    rst(6,8) = rst(8,6)
    rst(8,8) = rst(2,2)
    rst(12,8) = -rst(6,2)
    rst(3,9) = rst(9,3)
    rst(5,9) = rst(9,5)
    rst(9,9) = rst(3,3)
    rst(11,9) = -rst(5,3)
    rst(4,10) = rst(10,4)
    rst(10,10) = rst(4,4)
    rst(3,11) = rst(11,3)
    rst(5,11) = rst(11,5)
    rst(9,11) = rst(11,9)
    rst(11,11) = rst(5,5)
    rst(2,12) = rst(12,2)
    rst(6,12) = rst(12,6)
    rst(8,12) = rst(12,8)
    rst(12,12) = rst(6,6)

    ! Apply the transformation
    rst = matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
pure function b3d_mass_matrix(this, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(beam_element_3d), intent(in) :: this
        !! The beam_element_3d object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64) :: rho, L
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Initialization
    rho = this%material%density
    L = this%length()

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Compute the mass matrix
    allocate(rst(12, 12), source = 0.0d0)
    rst(1,1) = L * rho / 3.0d0
    rst(7,1) = L * rho / 6.0d0
    rst(2,2) = 1.3d1 * L * rho / 3.5d1
    rst(6,2) = 1.1d1 * rho * L**2 / 2.1d2
    rst(8,2) = 9.0d0 * L * rho / 7.0d1
    rst(12,2) = -1.3d1 * rho * L**2 / 4.2d2
    rst(3,3) = 1.3d0 * rho * L / 3.5d1
    rst(5,3) = -1.1d0 * rho * L**2 / 2.1d2
    rst(9,3) = 9.0d0 * rho * L / 7.0d1
    rst(11,3) = 1.3d1 * rho * L**2 / 4.2d2
    rst(4,4) = rho * L / 3.0d0
    rst(10,4) = rho * L / 6.0d0
    rst(3,5) = rst(5,3)
    rst(5,5) = rho * L**3 / 1.05d2
    rst(9,5) = -1.3d1 * rho * L**2 / 4.2d2
    rst(11,5) = -rho * L**3 / 1.4d2
    rst(2,6) = rst(6,2)
    rst(6,6) = rho * L**3 / 1.05d2
    rst(8,6) = 1.3d1 * rho * L**2 / 4.2d2
    rst(12,6) = -rho * L**3 / 1.4d2
    rst(1,7) = rst(7,1)
    rst(7,7) = rst(1,1)
    rst(2,8) = rst(8,2)
    rst(6,8) = rst(8,6)
    rst(8,8) = rst(2,2)
    rst(12,8) = -1.1d1 * rho * L**2 / 2.1d2
    rst(3,9) = rst(9,3)
    rst(5,9) = rst(9,5)
    rst(9,9) = rst(3,3)
    rst(11,9) = 1.1d1 * rho * L**2 / 2.1d2
    rst(4,10) = rst(10,4)
    rst(10,10) = rst(4,4)
    rst(3,11) = rst(11,3)
    rst(9,11) = rst(11,9)
    rst(11,11) = rst(5,5)
    rst(2,12) = rst(12,2)
    rst(6,12) = rst(12,6)
    rst(8,12) = rst(12,8)
    rst(12,12) = rst(6,6)

    ! Apply the transformation
    rst = matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
end module