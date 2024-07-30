! Shape Functions:
! 2D Line: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp1.html#shp2dlinerdof
! 3D Line: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp2.html#shp3d2node

module dynamics_structural
    use iso_fortran_env
    use linalg, only : csr_matrix, create_empty_csr_matrix
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
    type :: node
        !! Defines a node.
        integer(int32) :: index
            !! The global index of the node.
        real(real64) :: x
            !! The x-coordinate of the node.
        real(real64) :: y
            !! The y-coordinate of the node.
        real(real64) :: z
            !! The z-coordinate of the node.
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
        !! Defines a two-dimensional beam-element.
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
pure function create_connectivity_matrix(gdof, e, nodes) result(rst)
    !! Creates a connectivity matrix for the element.
    integer(int32), intent(in) :: gdof
        !! The number of global degrees of freedom.
    class(element), intent(in) :: e
        !! The element.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    type(csr_matrix) :: rst
        !! The resulting connectivity matrix.

    ! Local Variables
    integer(int32) :: i, j, col, nnodes, nnz, row
    integer(int32), allocatable, dimension(:) :: rows, cols
    real(real64), allocatable, dimension(:) :: vals

    ! Initialization
    nnodes = e%get_node_count()
    nnz = e%get_dof_per_node() * nnodes
    allocate(rows(nnz), cols(nnz), source = 0)
    allocate(vals(nnz), source = 1.0d0)

    ! Process
    row = 0
    do j = 1, nnodes
        col = find_global_dof(e%get_node(j), nodes)
        do i = 1, e%get_dof_per_node()
            row = row + 1
            rows(i) = row
            cols(i) = col
            col = col + 1
        end do
    end do

    ! TO DO: Build the matrix
end function

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
    rst(5,5) = st
    rst(6,6) = 1.0d0
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
end module