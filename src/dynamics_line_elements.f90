! Shape Functions:
! 2D Line: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp1.html#shp2dlinerdof
! 3D Line: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_shp2.html#shp3d2node

module dynamics_line_elements
    use iso_fortran_env
    use dynamics_structural
    use dynamics_geometry
    use dynamics_rotation
    use dynamics_error_handling
    implicit none
    private
    public :: beam_element_2d
    public :: beam_element_3d

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

    interface beam_element_2d
        module procedure :: b2d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(line_element) :: beam_element_3d
        !! Defines a three-dimensional Bernoulli-Euler beam element.
        real(real64) :: Ixx
            !! The beam moment of inertia about the element x-axis.
        real(real64) :: Iyy
            !! The beam moment of inertia about the element y-axis.
        real(real64) :: Izz
            !! The beam moment of inertia about the element z-axis.
        real(real64) :: Iyz
            !! The cross-sectional product of inertia.
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

    interface beam_element_3d
        module procedure :: b3d_init
    end interface

contains


! ******************************************************************************
! BEAM_2D ROUTINES
! ------------------------------------------------------------------------------
pure function b2d_init(mat, area, moi, nd1, nd2) result(rst)
    !! Initializes a new [[beam_element_2d]].
    class(material), intent(in) :: mat
        !! The material.
    real(real64), intent(in) :: area
        !! The cross-sectional area.
    real(real64), intent(in) :: moi
        !! The moment-of-inertia.
    class(node), intent(in) :: nd1
        !! Node 1 of the element.
    class(node), intent(in) :: nd2
        !! Node 2 of the element.
    type(beam_element_2d) :: rst
        !! The new [[beam_element_2d]].
    rst%material = mat
    rst%area = area
    rst%moment_of_inertia = moi
    rst%node_1 = nd1
    rst%node_2 = nd2
end function

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
    real(real64) :: l, sv
    
    ! Initialization
    allocate(rst(2, 6), source = 0.0d0)

    ! Process - use exact analytical formulas derived from the Hermite
    ! shape functions for a 2D Euler-Bernoulli beam element.
    ! For natural coordinate s in [-1, 1] and physical length l:
    !   Axial strain:  epsilon = du/dx
    !   Curvature:     kappa   = d2v/dx2
    l = this%length()
    sv = s(1)
    rst(1,1) = -1.0d0 / l
    rst(2,2) = 6.0d0 * sv / l**2
    rst(2,3) = (3.0d0 * sv - 1.0d0) / l
    rst(1,4) = 1.0d0 / l
    rst(2,5) = -6.0d0 * sv / l**2
    rst(2,6) = (3.0d0 * sv + 1.0d0) / l
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
! BEAM_3D ROUTINES
! ------------------------------------------------------------------------------
pure function b3d_init(mat, area, ixx, iyy, izz, iyz, nd1, nd2, orient) result(rst)
    !! Constructs a new [[beam_element_3d]] object.
    class(material), intent(in) :: mat
        !! The material.
    real(real64), intent(in) :: area
        !! The cross-sectional area.
    real(real64), intent(in) :: ixx
        !! The moment of inertia about the element's x-axis (torsional).
    real(real64), intent(in) :: iyy
        !! The moment of inertia about the element's y-axis.
    real(real64), intent(in) :: izz
        !! The moment of inertia about the element's z-axis.
    real(real64), intent(in) :: iyz
        !! The cross-sectional product of inertia.
    class(node), intent(in) :: nd1
        !! The first node.
    class(node), intent(in) :: nd2
        !! The second node.
    class(point), intent(in) :: orient
        !! An orientation node that locates the direction of the element
        !! z-axis.  The orientation point is measured relative to the first
        !! node in the element.  Specifically, the element z-axis is assumed
        !! to be defined by the location of this point relative to the
        !! location of node 1.
    type(beam_element_3d) :: rst
        !! The new [[beam_element_3d]] object.

    rst%material = mat
    rst%area = area
    rst%Ixx = ixx
    rst%Iyy = iyy
    rst%Izz = izz
    rst%Iyz = iyz
    rst%node_1 = nd1
    rst%node_2 = nd2
    rst%orientation_point = orient
end function

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
    rst(2,3) = this%Iyz * this%material%modulus
    rst(3,3) = this%Iyy * this%material%modulus
    rst(3,2) = rst(2,3)
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
    real(real64) :: A, E, Iyy, Izz, Iyz, Jxx, G, L
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Initialization
    A = this%area
    Jxx = this%Ixx
    Iyy = this%Iyy
    Izz = this%Izz
    Iyz = this%Iyz
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
    rst(2,3) = 1.2d1 * E * Iyz / L**3
    rst(5,2) = -6.0d0 * E * Iyz / L**2
    rst(8,3) = -rst(2,3)
    rst(11,2) = 6.0d0 * E * Iyz / L**2
    rst(6,3) = 6.0d0 * E * Iyz / L**2
    rst(8,5) = 6.0d0 * E * Iyz / L**2
    rst(12,3) = -6.0d0 * E * Iyz / L**2
    rst(5,9) = 6.0d0 * E * Iyz / L**2
    rst(11,9) = -6.0d0 * E * Iyz / L**2
    rst(6,5) = -4.0d0 * E * Iyz / L
    rst(12,5) = -2.0d0 * E * Iyz / L
    rst(6,11) = -2.0d0 * E * Iyz / L
    rst(12,11) = -4.0d0 * E * Iyz / L
    rst(4,4) = G * Jxx / L
    rst(10,4) = -rst(4,4)
    rst(3,5) = rst(5,3)
    rst(2,5) = rst(5,2)
    rst(2,9) = -rst(2,3)
    rst(2,11) = rst(11,2)
    rst(3,6) = rst(6,3)
    rst(3,8) = rst(8,3)
    rst(3,12) = rst(12,3)
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
    rst(5,6) = rst(6,5)
    rst(5,8) = rst(8,5)
    rst(5,12) = rst(12,5)
    rst(6,8) = rst(8,6)
    rst(8,8) = rst(2,2)
    rst(12,8) = -rst(6,2)
    rst(3,9) = rst(9,3)
    rst(6,9) = -6.0d0 * E * Iyz / L**2
    rst(8,9) = -rst(2,9)
    rst(9,2) = rst(2,9)
    rst(9,6) = rst(6,9)
    rst(9,8) = rst(8,9)
    rst(5,9) = rst(9,5)
    rst(9,9) = rst(3,3)
    rst(11,9) = -rst(5,3)
    rst(4,10) = rst(10,4)
    rst(10,10) = rst(4,4)
    rst(3,11) = rst(11,3)
    rst(8,11) = -6.0d0 * E * Iyz / L**2
    rst(9,11) = rst(11,9)
    rst(11,6) = rst(6,11)
    rst(11,8) = rst(8,11)
    rst(5,11) = rst(11,5)
    rst(9,11) = rst(11,9)
    rst(11,11) = rst(5,5)
    rst(2,12) = rst(12,2)
    rst(9,12) = -6.0d0 * E * Iyz / L**2
    rst(12,9) = rst(9,12)
    rst(11,12) = rst(12,11)
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