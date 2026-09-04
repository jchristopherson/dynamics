module dynamics_parallel_linkage
    !! Provides types supporting the kinematic analysis of closed-loop, or
    !! parallel, linkages.  The topology of the mechanism is described by a
    !! graph whose vertices are the links and whose edges are the joints.  A
    !! spanning tree of that graph supplies the transformation path to each
    !! link, and the edges excluded from the tree define the loop-closure
    !! constraints that must be satisfied by any valid configuration.
    use iso_fortran_env
    use dynamics_error_handling
    use dynamics_kinematics
    use dynamics_graph
    use dynamics_joints
    use dynamics_linkage, only : link
    use linalg, only : identity
    implicit none
    private
    public :: link_container
    public :: kinematic_mechanism
    public :: parallel_linkage
    public :: planar_linkage

    real(real64), parameter :: fd_step = 1.0d-6
        ! The finite difference step size used when computing Jacobians.

    type link_container
        !! A container allowing a collection of differing link types to be
        !! stored within a single array.
        class(link), allocatable, public :: item
            !! The link object.
    end type

    type, abstract :: kinematic_mechanism
        !! Defines the behavior common to all multi-link mechanisms, including
        !! those containing closed kinematic loops.  The joint variables of the
        !! mechanism are collected into a single array; the variables belonging
        !! to each joint appear in the order in which the joints were supplied.
        type(link_container), allocatable, private, dimension(:) :: m_links
        type(joint), allocatable, private, dimension(:) :: m_joints
        real(real64), allocatable, private, dimension(:,:,:) :: m_parentFrames
        real(real64), allocatable, private, dimension(:,:,:) :: m_childFrames
        type(graph), private :: m_graph
        type(spanning_tree), private :: m_tree
        type(graph_loop), allocatable, private, dimension(:) :: m_loops
        integer(int32), allocatable, private, dimension(:) :: m_qIndex
        integer(int32), allocatable, private, dimension(:) :: m_actuated
        integer(int32), allocatable, private, dimension(:) :: m_passive
        real(real64), allocatable, private, dimension(:) :: m_q
        integer(int32), private :: m_base = 1
        integer(int32), private :: m_effector = 1
        real(real64), private :: m_tool(4, 4)
    contains
        procedure(mechanism_dimension), deferred, public :: get_space_dimension
        procedure(mechanism_pose_error), deferred, public :: pose_error
        procedure, public :: get_link_count => km_get_link_count
        procedure, public :: get_link => km_get_link
        procedure, public :: get_joint_count => km_get_joint_count
        procedure, public :: get_joint => km_get_joint
        procedure, public :: get_variable_count => km_get_variable_count
        procedure, public :: get_joint_variable_index => km_get_variable_index
        procedure, public :: get_loop_count => km_get_loop_count
        procedure, public :: get_constraint_count => km_get_constraint_count
        procedure, public :: get_degrees_of_freedom => km_get_dof
        procedure, public :: get_actuated_variable_count => km_get_nactuated
        procedure, public :: get_actuated_indices => km_get_actuated_indices
        procedure, public :: get_actuated_variables => km_get_actuated_variables
        procedure, public :: get_base_link => km_get_base
        procedure, public :: get_end_effector_link => km_get_effector
        procedure, public :: get_tool_frame => km_get_tool
        procedure, public :: set_tool_frame => km_set_tool
        procedure, public :: get_configuration => km_get_configuration
        procedure, public :: set_configuration => km_set_configuration
        procedure, public :: joint_transform => km_joint_transform
        procedure, public :: body_transform => km_body_transform
        procedure, public :: end_effector_transform => km_effector_transform
        procedure, public :: constraints => km_constraints
        procedure, public :: constraint_jacobian => km_constraint_jacobian
        procedure, public :: solve_configuration => km_solve_configuration
        procedure, public :: forward_kinematics => km_forward_kinematics
        procedure, public :: jacobian => km_jacobian
        procedure, public :: inverse_kinematics => km_inverse_kinematics
    end type

    interface
        pure function mechanism_dimension(this) result(rst)
            !! Gets the dimension of the space in which the mechanism operates.
            use iso_fortran_env, only : int32
            import kinematic_mechanism
            class(kinematic_mechanism), intent(in) :: this
                !! The mechanism object.
            integer(int32) :: rst
                !! The number of independent components describing the
                !! displacement of one body relative to another.
        end function

        pure function mechanism_pose_error(this, x) result(rst)
            !! Reduces a transformation matrix describing the deviation between
            !! two coordinate frames into a vector of independent errors.  The
            !! resulting vector is zero when the supplied matrix is the identity
            !! matrix.
            use iso_fortran_env, only : real64
            import kinematic_mechanism
            class(kinematic_mechanism), intent(in) :: this
                !! The mechanism object.
            real(real64), intent(in) :: x(4, 4)
                !! The 4-by-4 error transformation matrix.
            real(real64), allocatable, dimension(:) :: rst
                !! The resulting error vector.
        end function
    end interface

    type, extends(kinematic_mechanism) :: parallel_linkage
        !! Defines a spatial, closed-loop linkage.  Each independent loop within
        !! the mechanism contributes six loop-closure constraint equations.
    contains
        procedure, public :: get_space_dimension => pl_get_space_dimension
        procedure, public :: pose_error => pl_pose_error
    end type

    interface parallel_linkage
        module procedure :: pl_init_containers
        module procedure :: pl_init_array
    end interface

    type, extends(kinematic_mechanism) :: planar_linkage
        !! Defines a planar, closed-loop linkage operating within the x-y plane.
        !! Each independent loop within the mechanism contributes three
        !! loop-closure constraint equations, avoiding the redundant equations
        !! that arise when a planar mechanism is analyzed in three dimensions.
        !!
        !! All joints must be usable within a planar mechanism, and all link
        !! geometry must be arranged such that motion remains within the x-y
        !! plane.
    contains
        procedure, public :: get_space_dimension => pln_get_space_dimension
        procedure, public :: pose_error => pln_pose_error
        procedure, public :: end_effector_pose => pln_end_effector_pose
    end type

    interface planar_linkage
        module procedure :: pln_init_containers
        module procedure :: pln_init_array
    end interface

    type mechanism_solver_data
        ! An internal container used to pass data to the nonlinear solver.
        class(kinematic_mechanism), pointer :: mechanism => null()
            ! The mechanism being solved.
        real(real64) :: target(4, 4) = 0.0d0
            ! The 4-by-4 end-effector target transformation matrix.
        real(real64), allocatable, dimension(:) :: actuated_targets
            ! The target values of the actuated joint variables.
        logical :: solve_for_pose = .false.
            ! True to drive the end-effector to the target transformation
            ! matrix; false to drive the actuated variables to their targets.
    end type

contains
! ******************************************************************************
! CONSTRUCTION
! ------------------------------------------------------------------------------
    subroutine initialize_mechanism(this, lnks, jnts, base, effector, tool)
        !! Establishes the topology of the mechanism and validates its
        !! description.
        class(kinematic_mechanism), intent(inout), target :: this
            !! The mechanism object.
        type(link_container), intent(in), dimension(:) :: lnks
            !! The collection of links forming the mechanism.
        type(joint), intent(in), dimension(:) :: jnts
            !! The collection of joints connecting the links.
        integer(int32), intent(in), optional :: base
            !! The index of the link that is fixed to ground.  If not supplied,
            !! the first link is used.
        integer(int32), intent(in), optional :: effector
            !! The index of the link carrying the end-effector.  If not
            !! supplied, the last link is used.
        real(real64), intent(in), optional :: tool(4, 4)
            !! The transformation matrix relating the end-effector coordinate
            !! frame to the body frame of the end-effector link.  If not
            !! supplied, an identity matrix is used.

        ! Local Variables
        integer(int32) :: i, j, k, m, n, nl, nj, ndof, nact, npassive
        class(link), pointer :: lnk

        ! Initialization
        nl = size(lnks)
        nj = size(jnts)
        if (nl < 1 .or. nj < 1) error stop DYN_INVALID_INPUT_ERROR

        this%m_base = 1
        this%m_effector = nl
        this%m_tool = identity(4)
        if (present(base)) this%m_base = base
        if (present(effector)) this%m_effector = effector
        if (present(tool)) this%m_tool = tool
        if (this%m_base < 1 .or. this%m_base > nl) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if
        if (this%m_effector < 1 .or. this%m_effector > nl) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if

        ! Store the links
        allocate(this%m_links(nl))
        do i = 1, nl
            if (.not.allocated(lnks(i)%item)) error stop DYN_NULL_POINTER_ERROR
            allocate(this%m_links(i)%item, source = lnks(i)%item)
        end do

        ! Store the joints and resolve the attachment frames
        allocate(this%m_joints(nj), source = jnts)
        allocate(this%m_parentFrames(4, 4, nj), this%m_childFrames(4, 4, nj))
        do i = 1, nj
            if (jnts(i)%parent_link < 1 .or. jnts(i)%parent_link > nl) then
                error stop DYN_INDEX_OUT_OF_RANGE
            end if
            if (jnts(i)%child_link < 1 .or. jnts(i)%child_link > nl) then
                error stop DYN_INDEX_OUT_OF_RANGE
            end if
            if (jnts(i)%parent_link == jnts(i)%child_link) then
                error stop DYN_INVALID_INPUT_ERROR
            end if
            lnk => this%m_links(jnts(i)%parent_link)%item
            if (jnts(i)%parent_frame > lnk%get_joint_count()) then
                error stop DYN_INDEX_OUT_OF_RANGE
            end if
            this%m_parentFrames(:,:,i) = lnk%get_joint_frame(jnts(i)%parent_frame)
            lnk => this%m_links(jnts(i)%child_link)%item
            if (jnts(i)%child_frame > lnk%get_joint_count()) then
                error stop DYN_INDEX_OUT_OF_RANGE
            end if
            this%m_childFrames(:,:,i) = lnk%get_joint_frame(jnts(i)%child_frame)
        end do

        ! Build the topology
        call this%m_graph%initialize(nl, nj)
        do i = 1, nj
            call this%m_graph%add_edge(jnts(i)%parent_link, jnts(i)%child_link)
        end do
        this%m_tree = this%m_graph%build_spanning_tree(this%m_base)
        if (.not.this%m_tree%is_connected()) error stop DYN_INVALID_INPUT_ERROR
        this%m_loops = this%m_graph%find_independent_loops(this%m_tree)

        ! Map the joint variables
        allocate(this%m_qIndex(nj + 1))
        this%m_qIndex(1) = 1
        do i = 1, nj
            ndof = jnts(i)%get_dof()
            if (ndof < 0) error stop DYN_INVALID_INPUT_ERROR
            this%m_qIndex(i + 1) = this%m_qIndex(i) + ndof
        end do
        n = this%m_qIndex(nj + 1) - 1
        if (n < 1) error stop DYN_INVALID_INPUT_ERROR
        allocate(this%m_q(n), source = 0.0d0)

        ! Identify the actuated and passive variables
        nact = 0
        do i = 1, nj
            if (jnts(i)%actuated) nact = nact + jnts(i)%get_dof()
        end do
        npassive = n - nact
        allocate(this%m_actuated(nact), this%m_passive(npassive))
        j = 0
        k = 0
        do i = 1, nj
            do m = this%m_qIndex(i), this%m_qIndex(i + 1) - 1
                if (jnts(i)%actuated) then
                    j = j + 1
                    this%m_actuated(j) = m
                else
                    k = k + 1
                    this%m_passive(k) = m
                end if
            end do
        end do
    end subroutine

! ------------------------------------------------------------------------------
    function to_containers(lnks) result(rst)
        ! Converts an array of links into an array of link containers.
        class(link), intent(in), dimension(:) :: lnks
        type(link_container), allocatable, dimension(:) :: rst
        integer(int32) :: i
        allocate(rst(size(lnks)))
        do i = 1, size(lnks)
            allocate(rst(i)%item, source = lnks(i))
        end do
    end function

! ------------------------------------------------------------------------------
    function pl_init_containers(lnks, jnts, base, effector, tool) result(rst)
        !! Initializes a new parallel_linkage object.
        type(link_container), intent(in), dimension(:) :: lnks
            !! The collection of links forming the mechanism.
        type(joint), intent(in), dimension(:) :: jnts
            !! The collection of joints connecting the links.
        integer(int32), intent(in), optional :: base
            !! The index of the link that is fixed to ground.  If not supplied,
            !! the first link is used.
        integer(int32), intent(in), optional :: effector
            !! The index of the link carrying the end-effector.  If not
            !! supplied, the last link is used.
        real(real64), intent(in), optional :: tool(4, 4)
            !! The transformation matrix relating the end-effector coordinate
            !! frame to the body frame of the end-effector link.  If not
            !! supplied, an identity matrix is used.
        type(parallel_linkage) :: rst
            !! The resulting parallel_linkage object.

        call initialize_mechanism(rst, lnks, jnts, base, effector, tool)
    end function

! ------------------------------------------------------------------------------
    function pl_init_array(lnks, jnts, base, effector, tool) result(rst)
        !! Initializes a new parallel_linkage object.
        class(link), intent(in), dimension(:) :: lnks
            !! The collection of links forming the mechanism.
        type(joint), intent(in), dimension(:) :: jnts
            !! The collection of joints connecting the links.
        integer(int32), intent(in), optional :: base
            !! The index of the link that is fixed to ground.  If not supplied,
            !! the first link is used.
        integer(int32), intent(in), optional :: effector
            !! The index of the link carrying the end-effector.  If not
            !! supplied, the last link is used.
        real(real64), intent(in), optional :: tool(4, 4)
            !! The transformation matrix relating the end-effector coordinate
            !! frame to the body frame of the end-effector link.  If not
            !! supplied, an identity matrix is used.
        type(parallel_linkage) :: rst
            !! The resulting parallel_linkage object.

        call initialize_mechanism(rst, to_containers(lnks), jnts, base, &
            effector, tool)
    end function

! ------------------------------------------------------------------------------
    function pln_init_containers(lnks, jnts, base, effector, tool) result(rst)
        !! Initializes a new planar_linkage object.
        type(link_container), intent(in), dimension(:) :: lnks
            !! The collection of links forming the mechanism.
        type(joint), intent(in), dimension(:) :: jnts
            !! The collection of joints connecting the links.  Every joint must
            !! be usable within a planar mechanism.
        integer(int32), intent(in), optional :: base
            !! The index of the link that is fixed to ground.  If not supplied,
            !! the first link is used.
        integer(int32), intent(in), optional :: effector
            !! The index of the link carrying the end-effector.  If not
            !! supplied, the last link is used.
        real(real64), intent(in), optional :: tool(4, 4)
            !! The transformation matrix relating the end-effector coordinate
            !! frame to the body frame of the end-effector link.  If not
            !! supplied, an identity matrix is used.
        type(planar_linkage) :: rst
            !! The resulting planar_linkage object.

        integer(int32) :: i
        do i = 1, size(jnts)
            if (.not.jnts(i)%is_planar_compatible()) then
                error stop DYN_INVALID_INPUT_ERROR
            end if
        end do
        call initialize_mechanism(rst, lnks, jnts, base, effector, tool)
    end function

! ------------------------------------------------------------------------------
    function pln_init_array(lnks, jnts, base, effector, tool) result(rst)
        !! Initializes a new planar_linkage object.
        class(link), intent(in), dimension(:) :: lnks
            !! The collection of links forming the mechanism.
        type(joint), intent(in), dimension(:) :: jnts
            !! The collection of joints connecting the links.  Every joint must
            !! be usable within a planar mechanism.
        integer(int32), intent(in), optional :: base
            !! The index of the link that is fixed to ground.  If not supplied,
            !! the first link is used.
        integer(int32), intent(in), optional :: effector
            !! The index of the link carrying the end-effector.  If not
            !! supplied, the last link is used.
        real(real64), intent(in), optional :: tool(4, 4)
            !! The transformation matrix relating the end-effector coordinate
            !! frame to the body frame of the end-effector link.  If not
            !! supplied, an identity matrix is used.
        type(planar_linkage) :: rst
            !! The resulting planar_linkage object.

        rst = pln_init_containers(to_containers(lnks), jnts, base, effector, &
            tool)
    end function

! ******************************************************************************
! PARALLEL_LINKAGE MEMBERS
! ------------------------------------------------------------------------------
    pure function pl_get_space_dimension(this) result(rst)
        !! Gets the dimension of the space in which the mechanism operates.
        class(parallel_linkage), intent(in) :: this
            !! The parallel_linkage object.
        integer(int32) :: rst
            !! The value six, corresponding to three translations and three
            !! rotations.

        rst = 6
    end function

! ------------------------------------------------------------------------------
    pure function pl_pose_error(this, x) result(rst)
        !! Reduces a transformation matrix describing the deviation between two
        !! coordinate frames into a six-element error vector.  The first three
        !! terms are the translational error, and the remaining three terms are
        !! an angle-axis approximation of the rotational error.
        class(parallel_linkage), intent(in) :: this
            !! The parallel_linkage object.
        real(real64), intent(in) :: x(4, 4)
            !! The 4-by-4 error transformation matrix.
        real(real64), allocatable, dimension(:) :: rst
            !! The resulting six-element error vector.

        allocate(rst(6))
        rst(1:3) = x(1:3,4)
        rst(4) = 0.5d0 * (x(3,2) - x(2,3))
        rst(5) = 0.5d0 * (x(1,3) - x(3,1))
        rst(6) = 0.5d0 * (x(2,1) - x(1,2))
    end function

! ******************************************************************************
! PLANAR_LINKAGE MEMBERS
! ------------------------------------------------------------------------------
    pure function pln_get_space_dimension(this) result(rst)
        !! Gets the dimension of the space in which the mechanism operates.
        class(planar_linkage), intent(in) :: this
            !! The planar_linkage object.
        integer(int32) :: rst
            !! The value three, corresponding to two translations and one
            !! rotation.

        rst = 3
    end function

! ------------------------------------------------------------------------------
    pure function pln_pose_error(this, x) result(rst)
        !! Reduces a transformation matrix describing the deviation between two
        !! coordinate frames into a three-element error vector containing the
        !! in-plane translational errors and the rotation about the plane
        !! normal.
        class(planar_linkage), intent(in) :: this
            !! The planar_linkage object.
        real(real64), intent(in) :: x(4, 4)
            !! The 4-by-4 error transformation matrix.
        real(real64), allocatable, dimension(:) :: rst
            !! The resulting three-element error vector.

        allocate(rst(3))
        rst(1) = x(1,4)
        rst(2) = x(2,4)
        rst(3) = atan2(x(2,1), x(1,1))
    end function

! ------------------------------------------------------------------------------
    function pln_end_effector_pose(this, q) result(rst)
        !! Computes the planar pose of the end-effector.
        class(planar_linkage), intent(inout), target :: this
            !! The planar_linkage object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing the actuated joint variables.
        real(real64) :: rst(3)
            !! A three-element array containing the x and y coordinates of the
            !! end-effector along with its orientation angle, in radians.

        real(real64) :: T(4, 4)
        T = this%forward_kinematics(q)
        rst(1) = T(1,4)
        rst(2) = T(2,4)
        rst(3) = atan2(T(2,1), T(1,1))
    end function

! ******************************************************************************
! KINEMATIC_MECHANISM MEMBERS
! ------------------------------------------------------------------------------
    pure function km_get_link_count(this) result(rst)
        !! Gets the number of links in the mechanism.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The link count.

        if (allocated(this%m_links)) then
            rst = size(this%m_links)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    function km_get_link(this, i) result(rst)
        !! Gets a pointer to the requested link object.
        class(kinematic_mechanism), intent(in), target :: this
            !! The mechanism object.
        integer(int32), intent(in) :: i
            !! The index of the link to retrieve (1 = first link).
        class(link), pointer :: rst
            !! A pointer to the requested link.

        if (i < 1 .or. i > this%get_link_count()) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if
        rst => this%m_links(i)%item
    end function

! ------------------------------------------------------------------------------
    pure function km_get_joint_count(this) result(rst)
        !! Gets the number of joints in the mechanism.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The joint count.

        if (allocated(this%m_joints)) then
            rst = size(this%m_joints)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    function km_get_joint(this, i) result(rst)
        !! Gets the requested joint object.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32), intent(in) :: i
            !! The index of the joint to retrieve (1 = first joint).
        type(joint) :: rst
            !! The requested joint.

        if (i < 1 .or. i > this%get_joint_count()) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if
        rst = this%m_joints(i)
    end function

! ------------------------------------------------------------------------------
    pure function km_get_variable_count(this) result(rst)
        !! Gets the total number of joint variables describing the configuration
        !! of the mechanism.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The number of joint variables.

        if (allocated(this%m_q)) then
            rst = size(this%m_q)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    function km_get_variable_index(this, i) result(rst)
        !! Gets the index, within the mechanism's array of joint variables, of
        !! the first variable belonging to the requested joint.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32), intent(in) :: i
            !! The index of the joint of interest (1 = first joint).
        integer(int32) :: rst
            !! The index of the joint's first variable.

        if (i < 1 .or. i > this%get_joint_count()) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if
        rst = this%m_qIndex(i)
    end function

! ------------------------------------------------------------------------------
    pure function km_get_loop_count(this) result(rst)
        !! Gets the number of independent kinematic loops in the mechanism.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The number of independent loops.

        if (allocated(this%m_loops)) then
            rst = size(this%m_loops)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    pure function km_get_constraint_count(this) result(rst)
        !! Gets the number of loop-closure constraint equations imposed upon the
        !! mechanism.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The number of constraint equations.

        rst = this%get_loop_count() * this%get_space_dimension()
    end function

! ------------------------------------------------------------------------------
    pure function km_get_dof(this) result(rst)
        !! Gets the number of degrees of freedom, or mobility, of the mechanism.
        !! The value is the number of joint variables less the number of
        !! loop-closure constraint equations.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The number of degrees of freedom.

        rst = this%get_variable_count() - this%get_constraint_count()
    end function

! ------------------------------------------------------------------------------
    pure function km_get_nactuated(this) result(rst)
        !! Gets the number of actuated joint variables.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The number of actuated joint variables.

        if (allocated(this%m_actuated)) then
            rst = size(this%m_actuated)
        else
            rst = 0
        end if
    end function

! ------------------------------------------------------------------------------
    function km_get_actuated_indices(this) result(rst)
        !! Gets the indices, within the mechanism's array of joint variables, of
        !! the actuated joint variables.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32), allocatable, dimension(:) :: rst
            !! An array containing the indices of the actuated variables.

        rst = this%m_actuated
    end function

! ------------------------------------------------------------------------------
    function km_get_actuated_variables(this, q) result(rst)
        !! Extracts the actuated joint variables from a full set of joint
        !! variables.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.
        real(real64), allocatable, dimension(:) :: rst
            !! An array containing the actuated joint variables.

        if (size(q) /= this%get_variable_count()) error stop DYN_ARRAY_SIZE_ERROR
        rst = q(this%m_actuated)
    end function

! ------------------------------------------------------------------------------
    pure function km_get_base(this) result(rst)
        !! Gets the index of the link that is fixed to ground.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The index of the base link.

        rst = this%m_base
    end function

! ------------------------------------------------------------------------------
    pure function km_get_effector(this) result(rst)
        !! Gets the index of the link carrying the end-effector.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32) :: rst
            !! The index of the end-effector link.

        rst = this%m_effector
    end function

! ------------------------------------------------------------------------------
    pure function km_get_tool(this) result(rst)
        !! Gets the transformation matrix relating the end-effector coordinate
        !! frame to the body frame of the end-effector link.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        real(real64) :: rst(4, 4)
            !! The 4-by-4 transformation matrix.

        rst = this%m_tool
    end function

! ------------------------------------------------------------------------------
    subroutine km_set_tool(this, x)
        !! Sets the transformation matrix relating the end-effector coordinate
        !! frame to the body frame of the end-effector link.
        class(kinematic_mechanism), intent(inout) :: this
            !! The mechanism object.
        real(real64), intent(in) :: x(4, 4)
            !! The 4-by-4 transformation matrix.

        this%m_tool = x
    end subroutine

! ------------------------------------------------------------------------------
    function km_get_configuration(this) result(rst)
        !! Gets the most recently computed set of joint variables.  This value
        !! is used as the starting estimate for subsequent solutions.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        real(real64), allocatable, dimension(:) :: rst
            !! An array containing all of the mechanism's joint variables.

        rst = this%m_q
    end function

! ------------------------------------------------------------------------------
    subroutine km_set_configuration(this, q)
        !! Sets the set of joint variables used as the starting estimate for
        !! subsequent solutions.  A reasonable estimate is important as a
        !! closed-loop mechanism typically admits multiple assembly modes.
        class(kinematic_mechanism), intent(inout) :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.

        if (size(q) /= this%get_variable_count()) error stop DYN_ARRAY_SIZE_ERROR
        this%m_q = q
    end subroutine

! ------------------------------------------------------------------------------
    function km_joint_transform(this, i, q) result(rst)
        !! Computes the transformation matrix relating the body coordinate frame
        !! of a joint's child link to the body coordinate frame of the joint's
        !! parent link.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32), intent(in) :: i
            !! The index of the joint of interest (1 = first joint).
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Local Variables
        integer(int32) :: i1, i2
        real(real64) :: M(4, 4)

        ! Input Checking
        if (i < 1 .or. i > this%get_joint_count()) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if
        if (size(q) /= this%get_variable_count()) error stop DYN_ARRAY_SIZE_ERROR

        ! Process
        i1 = this%m_qIndex(i)
        i2 = this%m_qIndex(i + 1) - 1
        M = this%m_joints(i)%motion_transform(q(i1:i2))
        rst = matmul(this%m_parentFrames(:,:,i), &
            matmul(M, transform_inverse(this%m_childFrames(:,:,i))))
    end function

! ------------------------------------------------------------------------------
    function km_body_transform(this, i, q) result(rst)
        !! Computes the transformation matrix relating the body coordinate frame
        !! of the requested link to the coordinate frame of the base link.  The
        !! transformation is accumulated along the path through the mechanism's
        !! spanning tree.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        integer(int32), intent(in) :: i
            !! The index of the link of interest (1 = first link).
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Local Variables
        integer(int32) :: k
        real(real64) :: T(4, 4)
        type(graph_path) :: path

        ! Input Checking
        if (i < 1 .or. i > this%get_link_count()) then
            error stop DYN_INDEX_OUT_OF_RANGE
        end if

        ! Process
        path = this%m_tree%get_path(i)
        rst = identity(4)
        do k = 1, size(path%edges)
            T = this%joint_transform(path%edges(k), q)
            if (path%forward(k)) then
                rst = matmul(rst, T)
            else
                rst = matmul(rst, transform_inverse(T))
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function km_effector_transform(this, q) result(rst)
        !! Computes the transformation matrix relating the end-effector
        !! coordinate frame to the coordinate frame of the base link.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        rst = matmul(this%body_transform(this%m_effector, q), this%m_tool)
    end function

! ------------------------------------------------------------------------------
    function km_constraints(this, q) result(rst)
        !! Evaluates the loop-closure constraint equations for the mechanism.  A
        !! valid configuration drives every equation to zero.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.
        real(real64), allocatable, dimension(:) :: rst
            !! An array containing the residual of each constraint equation.

        ! Local Variables
        integer(int32) :: i, i1, i2, nd
        real(real64) :: Tu(4, 4), Tv(4, 4), Te(4, 4), D(4, 4)

        ! Input Checking
        if (size(q) /= this%get_variable_count()) error stop DYN_ARRAY_SIZE_ERROR

        ! Process
        nd = this%get_space_dimension()
        allocate(rst(this%get_constraint_count()))
        do i = 1, this%get_loop_count()
            Tu = this%body_transform(this%m_loops(i)%vertex_1, q)
            Tv = this%body_transform(this%m_loops(i)%vertex_2, q)
            Te = this%joint_transform(this%m_loops(i)%cut_edge, q)
            D = matmul(transform_inverse(matmul(Tu, Te)), Tv)
            i1 = nd * (i - 1) + 1
            i2 = nd * i
            rst(i1:i2) = this%pose_error(D)
        end do
    end function

! ------------------------------------------------------------------------------
    function km_constraint_jacobian(this, q) result(rst)
        !! Computes the Jacobian matrix of the loop-closure constraint equations
        !! with respect to the mechanism's joint variables.  The derivatives are
        !! estimated by means of a central difference approximation.
        class(kinematic_mechanism), intent(in) :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing all of the mechanism's joint variables.
        real(real64), allocatable, dimension(:,:) :: rst
            !! The resulting C-by-N matrix where C is the number of constraint
            !! equations and N is the number of joint variables.

        ! Local Variables
        integer(int32) :: i, n, nc
        real(real64) :: h
        real(real64), allocatable, dimension(:) :: work, f1, f2

        ! Input Checking
        if (size(q) /= this%get_variable_count()) error stop DYN_ARRAY_SIZE_ERROR

        ! Process
        n = this%get_variable_count()
        nc = this%get_constraint_count()
        allocate(rst(nc, n))
        work = q
        do i = 1, n
            h = fd_step * max(abs(q(i)), 1.0d0)
            work(i) = q(i) + h
            f1 = this%constraints(work)
            work(i) = q(i) - h
            f2 = this%constraints(work)
            work(i) = q(i)
            rst(:,i) = (f1 - f2) / (2.0d0 * h)
        end do
    end function

! ------------------------------------------------------------------------------
    function km_solve_configuration(this, q, qo, ib) result(rst)
        !! Determines the complete set of joint variables satisfying both the
        !! loop-closure constraints and the requested actuated joint variables.
        class(kinematic_mechanism), intent(inout), target :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing the actuated joint variables.  The array must
            !! contain one value for each degree of freedom of the mechanism.
        real(real64), intent(in), optional, dimension(:) :: qo
            !! An optional array containing an initial estimate of all of the
            !! mechanism's joint variables.  If not supplied, the most recently
            !! computed configuration is used.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output that can be used to gather information on the
            !! solver.
        real(real64), allocatable, dimension(:) :: rst
            !! An array containing all of the mechanism's joint variables.

        ! Local Variables
        integer(int32) :: nc, na, nvar
        real(real64), allocatable, dimension(:) :: x
        type(mechanism_solver_data) :: obj

        ! Initialization
        nvar = this%get_variable_count()
        nc = this%get_constraint_count()
        na = this%get_actuated_variable_count()

        ! Input Checking
        if (na /= this%get_degrees_of_freedom()) error stop DYN_CONSTRAINT_ERROR
        if (size(q) /= na) error stop DYN_ARRAY_SIZE_ERROR
        if (present(qo)) then
            if (size(qo) /= nvar) error stop DYN_ARRAY_SIZE_ERROR
            x = qo
        else
            x = this%m_q
        end if

        ! Assemble and solve the system
        obj%mechanism => this
        obj%actuated_targets = q
        obj%solve_for_pose = .false.
        call solve_mechanism(x, nc + na, obj, ib)
        this%m_q = x
        rst = x
    end function

! ----------
    subroutine solve_mechanism(x, neqn, obj, ib)
        ! Solves the nonlinear system describing the configuration of the
        ! mechanism.  A Jacobian matrix is supplied to the solver as the default
        ! finite difference estimate is unreliable when a joint variable is
        ! near zero.
        real(real64), intent(inout), dimension(:) :: x
        integer(int32), intent(in) :: neqn
        type(mechanism_solver_data), intent(inout), target :: obj
        type(iteration_behavior), intent(out), optional :: ib

        real(real64), allocatable, dimension(:) :: resid
        type(vecfcn_helper) :: helper
        type(least_squares_solver) :: solver
        procedure(vecfcn), pointer :: fcn
        procedure(jacobianfcn), pointer :: jfcn

        if (neqn < size(x)) error stop DYN_CONSTRAINT_ERROR
        fcn => mechanism_equations
        jfcn => mechanism_equation_jacobian
        call helper%set_fcn(fcn, neqn, size(x))
        call helper%set_jacobian(jfcn)
        allocate(resid(neqn))
        call solver%solve(helper, x, resid, ib = ib, args = obj)
    end subroutine

! ----------
    subroutine mechanism_equations(x, f, args)
        ! The loop-closure equations along with either the actuation equations
        ! or the end-effector equations.
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(inout), optional :: args

        integer(int32) :: nc
        real(real64) :: T(4, 4)

        select type (args)
        class is (mechanism_solver_data)
            nc = args%mechanism%get_constraint_count()
            f(1:nc) = args%mechanism%constraints(x)
            if (args%solve_for_pose) then
                T = args%mechanism%end_effector_transform(x)
                f(nc+1:) = args%mechanism%pose_error( &
                    matmul(transform_inverse(T), args%target))
            else
                f(nc+1:) = x(args%mechanism%m_actuated) - args%actuated_targets
            end if
        end select
    end subroutine

! ----------
    subroutine mechanism_equation_jacobian(x, jac, args)
        ! The Jacobian matrix of the mechanism equations, estimated by means of
        ! a central difference approximation employing a step size that remains
        ! well conditioned as a joint variable approaches zero.
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(inout), optional :: args

        integer(int32) :: i, n, m
        real(real64) :: h
        real(real64), allocatable, dimension(:) :: work, f1, f2

        n = size(x)
        m = size(jac, 1)
        allocate(f1(m), f2(m))
        work = x
        do i = 1, n
            h = fd_step * max(abs(x(i)), 1.0d0)
            work(i) = x(i) + h
            call mechanism_equations(work, f1, args)
            work(i) = x(i) - h
            call mechanism_equations(work, f2, args)
            work(i) = x(i)
            jac(:,i) = (f1 - f2) / (2.0d0 * h)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    function km_forward_kinematics(this, q, qo, ib) result(rst)
        !! Computes the forward kinematics for the mechanism resulting in a
        !! transformation matrix relating the end-effector coordinate frame to
        !! the coordinate frame of the base link.
        !!
        !! Unlike a serial linkage, the forward kinematics of a closed-loop
        !! mechanism require the solution of the loop-closure constraints.  The
        !! solution obtained depends upon the assembly mode of the mechanism,
        !! and therefore upon the starting estimate.
        class(kinematic_mechanism), intent(inout), target :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing the actuated joint variables.  The array must
            !! contain one value for each degree of freedom of the mechanism.
        real(real64), intent(in), optional, dimension(:) :: qo
            !! An optional array containing an initial estimate of all of the
            !! mechanism's joint variables.  If not supplied, the most recently
            !! computed configuration is used.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output that can be used to gather information on the
            !! solver.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        real(real64), allocatable, dimension(:) :: x
        x = this%solve_configuration(q, qo, ib)
        rst = this%end_effector_transform(x)
    end function

! ------------------------------------------------------------------------------
    function km_jacobian(this, q, qo) result(rst)
        !! Constructs the Jacobian matrix relating the actuated joint velocities
        !! to the velocity of the end-effector, as expressed in the end-effector
        !! coordinate frame.
        !!
        !! The loop-closure constraints supply the relationship between the
        !! actuated and passive joint velocities.  Partitioning the constraint
        !! Jacobian into its actuated and passive terms yields
        !! $$ \dot{\vec{q}}_{p} = -C_{p}^{-1} C_{a} \dot{\vec{q}}_{a}, $$
        !! which is then combined with the end-effector Jacobian.  The partition
        !! becomes singular at a configuration in which the mechanism loses
        !! control of one or more of its degrees of freedom.
        class(kinematic_mechanism), intent(inout), target :: this
            !! The mechanism object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing the actuated joint variables.  The array must
            !! contain one value for each degree of freedom of the mechanism.
        real(real64), intent(in), optional, dimension(:) :: qo
            !! An optional array containing an initial estimate of all of the
            !! mechanism's joint variables.  If not supplied, the most recently
            !! computed configuration is used.
        real(real64), allocatable, dimension(:,:) :: rst
            !! The resulting D-by-M matrix where D is the dimension of the space
            !! in which the mechanism operates and M is the number of actuated
            !! joint variables.

        ! Local Variables
        integer(int32) :: nc, np
        real(real64), allocatable, dimension(:) :: x
        real(real64), allocatable, dimension(:,:) :: Cq, Cp, Ca, Je

        ! Initialization
        nc = this%get_constraint_count()
        np = size(this%m_passive)
        if (nc /= np) error stop DYN_CONSTRAINT_ERROR

        ! Establish the configuration of the mechanism
        x = this%solve_configuration(q, qo)

        ! Partition the constraint Jacobian and solve for the passive velocities
        Cq = this%constraint_jacobian(x)
        Cp = Cq(:,this%m_passive)
        Ca = -Cq(:,this%m_actuated)
        call solve_square_system(Cp, Ca)

        ! Combine with the end-effector Jacobian
        Je = end_effector_jacobian(this, x)
        rst = Je(:,this%m_actuated) + matmul(Je(:,this%m_passive), Ca)
    end function

! ----------
    subroutine solve_square_system(a, b)
        ! Solves the system A X = B by means of an LU factorization employing
        ! partial pivoting.  Both matrices are overwritten, with the solution
        ! being returned in B.
        real(real64), intent(inout), dimension(:,:) :: a
        real(real64), intent(inout), dimension(:,:) :: b

        integer(int32) :: i, j, k, n, ipiv
        real(real64) :: f, amax
        real(real64), allocatable, dimension(:) :: temp

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n) then
            error stop DYN_MATRIX_SIZE_ERROR
        end if
        allocate(temp(max(n, size(b, 2))))

        do k = 1, n
            ! Locate and apply the pivot
            ipiv = k
            amax = abs(a(k,k))
            do i = k + 1, n
                if (abs(a(i,k)) > amax) then
                    amax = abs(a(i,k))
                    ipiv = i
                end if
            end do
            ! A singular partition indicates the mechanism is at a singularity
            if (amax < epsilon(amax)) error stop DYN_CONSTRAINT_ERROR
            if (ipiv /= k) then
                temp(1:n) = a(k,:)
                a(k,:) = a(ipiv,:)
                a(ipiv,:) = temp(1:n)
                temp(1:size(b,2)) = b(k,:)
                b(k,:) = b(ipiv,:)
                b(ipiv,:) = temp(1:size(b,2))
            end if

            ! Eliminate
            do i = k + 1, n
                f = a(i,k) / a(k,k)
                a(i,k) = 0.0d0
                do j = k + 1, n
                    a(i,j) = a(i,j) - f * a(k,j)
                end do
                b(i,:) = b(i,:) - f * b(k,:)
            end do
        end do

        ! Back substitution
        do k = n, 1, -1
            do j = k + 1, n
                b(k,:) = b(k,:) - a(k,j) * b(j,:)
            end do
            b(k,:) = b(k,:) / a(k,k)
        end do
    end subroutine

! ----------
    function end_effector_jacobian(this, q) result(rst)
        ! Computes the derivative of the end-effector pose with respect to each
        ! of the mechanism's joint variables.
        class(kinematic_mechanism), intent(in) :: this
        real(real64), intent(in), dimension(:) :: q
        real(real64), allocatable, dimension(:,:) :: rst

        integer(int32) :: i, n, nd
        real(real64) :: h, T0i(4, 4), T1(4, 4), T2(4, 4)
        real(real64), allocatable, dimension(:) :: work

        n = this%get_variable_count()
        nd = this%get_space_dimension()
        allocate(rst(nd, n))
        T0i = transform_inverse(this%end_effector_transform(q))
        work = q
        do i = 1, n
            h = fd_step * max(abs(q(i)), 1.0d0)
            work(i) = q(i) + h
            T1 = this%end_effector_transform(work)
            work(i) = q(i) - h
            T2 = this%end_effector_transform(work)
            work(i) = q(i)
            rst(:,i) = (this%pose_error(matmul(T0i, T1)) - &
                this%pose_error(matmul(T0i, T2))) / (2.0d0 * h)
        end do
    end function

! ------------------------------------------------------------------------------
    function km_inverse_kinematics(this, trg, qo, ib) result(rst)
        !! Solves the inverse kinematics problem for the mechanism.  The
        !! loop-closure constraints are solved simultaneously with the
        !! end-effector constraints.
        class(kinematic_mechanism), intent(inout), target :: this
            !! The mechanism object.
        real(real64), intent(in) :: trg(4, 4)
            !! A transformation matrix relating the end-effector coordinate
            !! frame to the coordinate frame of the base link.  This
            !! transformation matrix defines the end-effector target for the
            !! solver.
        real(real64), intent(in), optional, dimension(:) :: qo
            !! An optional array containing an initial estimate of all of the
            !! mechanism's joint variables.  If not supplied, the most recently
            !! computed configuration is used.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output that can be used to gather information on the
            !! solver.
        real(real64), allocatable, dimension(:) :: rst
            !! An array containing the actuated joint variables that satisfy the
            !! constraints.

        ! Local Variables
        integer(int32) :: nc, nd, nvar
        real(real64), allocatable, dimension(:) :: x
        type(mechanism_solver_data) :: obj

        ! Initialization
        nvar = this%get_variable_count()
        nc = this%get_constraint_count()
        nd = this%get_space_dimension()

        ! Input Checking
        if (nd < this%get_degrees_of_freedom()) error stop DYN_CONSTRAINT_ERROR
        if (present(qo)) then
            if (size(qo) /= nvar) error stop DYN_ARRAY_SIZE_ERROR
            x = qo
        else
            x = this%m_q
        end if

        ! Assemble and solve the system
        obj%mechanism => this
        obj%target = trg
        obj%solve_for_pose = .true.
        call solve_mechanism(x, nc + nd, obj, ib)
        this%m_q = x
        rst = this%get_actuated_variables(x)
    end function

! ------------------------------------------------------------------------------
end module
