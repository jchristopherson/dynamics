module dynamics_linkage
    use iso_fortran_env
    use dynamics_rigid_bodies
    use dynamics_kinematics
    use dynamics_geometry
    use dynamics_helper
    use dynamics_quaternions
    use collections, only : list
    use ferror, only : errors
    implicit none
    private
    public :: binary_link
    public :: serial_linkage

    type, extends(rigid_body) :: binary_link
        !! Defines a link consisting of only two joints.  The coordinate system
        !! of this link is situated at the distal joint with it's z-axis 
        !! coincident with the axis of the joint.  The link utilizes a
        !! Denavit-Hartenberg convention in order to express its geometry.
        real(real64), public :: link_length
            !! The link length is the distance between the proximal and distal
            !! joint axes as measured along the link's x-axis.
        real(real64), public :: link_twist
            !! The link twist is the required rotation of the proximal joint 
            !! axis about the link's x-axis to become parallel to the distal
            !! joint's axis.
        real(real64), public :: link_offset
            !! The link offset is the fixed distance between the previous link's
            !! x-axis and the current link's x-axis as measured along the
            !! axis of the proximal joint.
        real(real64), public :: joint_angle
            !! The joint angle is the required rotation of the previous link's
            !! x-axis about the proximal joint's axis to become parallel to the
            !! current link's x-axis.
        real(real64), public :: joint_type
            !! The proximal joint type.  This value must be either 
            !! REVOLUTE_JOINT or PRISMATIC_JOINT.
    end type

    interface binary_link
        module procedure :: bl_init
    end interface

    type serial_linkage
        !! Defines a serial linkage.
        type(list), private :: m_links
    contains
        procedure, public :: get_link_count => sl_get_link_count
        procedure, public :: get_link => sl_get_link
        procedure, public :: forward_kinematics => sl_forward_kinematics
        procedure, public :: jacobian => sl_jacobian
        procedure, public :: inverse_kinematics => sl_inverse_kinematics_1
    end type

    interface serial_linkage
        module procedure :: sl_init
    end interface

    type serial_linkage_solver_data
        ! An internal type used as a container to pass data to the inverse
        ! kinematics solver.
        type(serial_linkage), pointer :: linkage
            ! The serial_linkage object.
        real(real64) :: target(4, 4)
            ! The 4-by-4 target transformation matrix.
    end type

contains
! ******************************************************************************
! BINARY_LINK MEMBERS
! ------------------------------------------------------------------------------
    pure function bl_init(jtype, length, twist, offset, angle, &
        mass, inertia, cg) result(rst)
        !! Initializes a new parallel_revolute_revolute_link instance.
        integer(int32), intent(in), optional :: jtype
            !! The proximal joint type.  This value must be either &
            !! REVOLUTE_JOINT or PRISMATIC_JOINT.  If incorrectly specified, 
            !! this parameter defaults to REVOLUTE_JOINT.
        real(real64), intent(in), optional :: length
            !! The link length.  If no value is specified, a value of 0 is used.
        real(real64), intent(in), optional :: twist
            !! The link twist angle.  If no value is specified, a value of 0
            !! is used.
        real(real64), intent(in), optional :: offset
            !! The link offset.  If no value is specified, a value of 0 is used.
        real(real64), intent(in), optional :: angle
            !! The joint angle offset.  If no value is specified, a value of 0
            !! is used.
        real(real64), intent(in), optional :: mass
            !! The mass of the link.  If no value is specified, a value of 1 is
            !! used.
        real(real64), intent(in), optional :: inertia(3, 3)
            !! The 3-by-3 inertia tensor.  If not specified, an identity matrix
            !! is used.
        real(real64), intent(in), optional :: cg(3)
            !! The x-y-z location of the CG relative to the distal coordinate
            !! frame, expressed in the link coordinate frame.  If not supplied,
            !! the CG is set to (0, 0, 0) such that it is located at the center
            !! of the distal joint.
        type(binary_link) :: rst
            !! The resulting binary_link object.

        ! Initialize the base object
        call initialize_rigid_body(rst, mass, inertia, cg)
        
        ! Additional initialization
        if (present(jtype)) then
            if (jtype /= REVOLUTE_JOINT .and. jtype /= PRISMATIC_JOINT) then
                rst%joint_type = REVOLUTE_JOINT
            else
                rst%joint_type = jtype
            end if
        else
            rst%joint_type = REVOLUTE_JOINT
        end if

        if (present(length)) then
            rst%link_length = length
        else
            rst%link_length = 0.0d0
        end if
        
        if (present(twist)) then
            rst%link_twist = twist
        else
            rst%link_twist = 0.0d0
        end if

        if (present(offset)) then
            rst%link_offset = offset
        else
            rst%link_offset = 0.0d0
        end if

        if (present(angle)) then
            rst%joint_angle = angle
        else
            rst%joint_angle = 0.0d0
        end if
    end function

! ******************************************************************************
! SERIAL_LINKAGE MEMBERS
! ------------------------------------------------------------------------------
    function sl_init(lnks) result(rst)
        !! Initializes a new serial_linkage object.
        class(binary_link), intent(in), dimension(:) :: lnks
            !! A collection of binary_link objects.  The collection starts with
            !! the first link and progresses to the end-effector in a 
            !! sequential manner.
        type(serial_linkage) :: rst
            !! The serial_linkage instance.

        ! Initialization
        integer(int32) :: i, n
        n = size(lnks)
        do i = 1, n
            call rst%m_links%push(lnks(i))
        end do
    end function

! ------------------------------------------------------------------------------
    pure function sl_get_link_count(this) result(rst)
        !! Gets the number of links in the linkage.
        class(serial_linkage), intent(in) :: this
            !! The serial_linkage object.
        integer(int32) :: rst
            !! The link count.

        rst = this%m_links%count()
    end function

! ------------------------------------------------------------------------------
    function sl_get_link(this, i) result(rst)
        !! Gets a pointer to the requested link object.
        class(serial_linkage), intent(in) :: this
            !! The serial_linkage object.
        integer(int32), intent(in) :: i
            !! The index of the link to retrieve (1 = first link).
        class(binary_link), pointer :: rst
            !! A pointer to the requested link.

        ! Process
        class(*), pointer :: ptr
        ptr => this%m_links%get(i)
        rst => null()
        select type (ptr)
        class is (binary_link)
            rst => ptr
        end select
    end function

! ------------------------------------------------------------------------------
    function sl_forward_kinematics(this, q, err) result(rst)
        use dynamics_error_handling, only : report_array_size_error
        !! Computes the forward kinematics for the linkage resulting in a 
        !! transformation matrix between world and end-effector coordinate
        !! frames.
        class(serial_linkage), intent(in) :: this
            !! The serial_linkage object.
        real(real64), intent(in), dimension(:) :: q
            !! The array of joint variables.  This array must be the same size
            !! as there are number of links in this linkage.
        class(errors), intent(inout), optional, target :: err
            !! An errors-based object providing error handling in the event the
            !! array size is incorrect.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: qt, qd, T(4, 4)
        class(binary_link), pointer :: lnk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = this%get_link_count()
        rst = identity_4()

        ! Input Checking
        if (size(q) /= n) then
            call report_array_size_error("sl_forward_kinematics", "q", n, &
                size(q), errmgr)
            return
        end if

        ! Process
        do i = 1, n
            ! Compute the link transformation matrix
            lnk => this%get_link(i)
            qd = lnk%link_offset
            qt = lnk%joint_angle
            if (lnk%joint_type == PRISMATIC_JOINT) then
                qd = qd + q(i)
            else
                qt = qt + q(i)
            end if
            T = dh_matrix(lnk%link_twist, lnk%link_length, qt, qd)

            ! Accumulate the transform
            rst = matmul(rst, T)
        end do
    end function

! ------------------------------------------------------------------------------
    function sl_jacobian(this, q, err) result(rst)
        use dynamics_error_handling, only : report_array_size_error
        !! Constructs the Jacobian matrix for the linkage.  The Jacobian matrix 
        !! relates the joint velocities \(\dot{\vec{q}}\) to the end-effector 
        !! velocity \(\dot{\vec{X}}\) by \(\dot{\vec{X}} = J \dot{\vec{q}}\).
        class(serial_linkage), intent(in) :: this
            !! The serial_linkage object.
        real(real64), intent(in), dimension(:) :: q
            !! The array of joint variables.  This array must be the same size
            !! as there are number of links in this linkage.
        class(errors), intent(inout), optional, target :: err
            !! An errors-based object providing error handling in the event the
            !! array size is incorrect.
        real(real64), allocatable, dimension(:,:) :: rst
            !! The resulting 6-by-N Jacobian matrix where N is the number of
            !! links in the linkage.

        ! Parameters
        real(real64), parameter :: zi_1(4) = [0.0d0, 0.0d0, 1.0d0, 0.0d0]

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: qd, qt
        integer(int32), allocatable, dimension(:) :: jtypes
        real(real64), allocatable, dimension(:) :: a, alpha, theta, d
        class(binary_link), pointer :: lnk
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = this%get_link_count()

        ! Error Checking
        if (size(q) /= n) then
            call report_array_size_error("sl_jacobian", "q", n, size(q), errmgr)
            return
        end if

        ! Format inputs for the Jacobian calculation
        allocate(a(n), alpha(n), theta(n), d(n), jtypes(n))
        do i = 1, n
            lnk => this%get_link(i)
            qd = lnk%link_offset
            qt = lnk%joint_angle
            if (lnk%joint_type == PRISMATIC_JOINT) then
                qd = qd + q(i)
            else
                qt = qt + q(i)
            end if
            a(i) = lnk%link_length
            alpha(i) = lnk%link_twist
            theta(i) = qt
            d(i) = qd
            jtypes(i) = lnk%joint_type
        end do

        ! Compute the Jacobian
        rst = dh_jacobian(alpha, a, theta, d, jtypes)
    end function

! ------------------------------------------------------------------------------
    function sl_inverse_kinematics_1(this, qo, trg, ib, err) result(rst)
        !! Solves the inverse kinematics problem for the linkage.
        class(serial_linkage), intent(in), target :: this
            !! The serial_linkage object.
        real(real64), intent(in), dimension(:) :: qo
            !! An M-element array containing an initial estimate of the M joint
            !! variables.
        real(real64), intent(in) :: trg(4, 4)
            !! A transformation matrix relating the end-effector coordinate
            !! frame and the world coordinate frame.  This transformation 
            !! matrix defines the end-effector target for the solver.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output that can be used to gather information on the
            !! solver.
        class(errors), intent(inout), optional, target :: err
            !! An optional error handling object used to retrieve any errors
            !! regarding the solver.
        real(real64), allocatable, dimension(:) :: rst
            !! An M-element array containing the computed joint variables that
            !! satisfy the constraints.

        ! Local Variables
        procedure(vecfcn), pointer :: vfcn
        procedure(jacobianfcn), pointer :: jfcn
        type(serial_linkage_solver_data) :: obj
        real(real64) :: constraints(6)

        ! Initialization
        vfcn => sl_vecfcn
        jfcn => sl_jacobianfcn
        obj%linkage => this
        obj%target = trg
        constraints(1:3) = trg(1:3,4)
        constraints(4:6) = 0.0d0

        ! Input Check
        ! TO DO: If there are more than 6 joint variables we'll need to error
        ! out for now until we develop a better solver

        ! Process
        rst = solve_inverse_kinematics(vfcn, qo, constraints, ib = ib, &
            args = obj, jfcn = jfcn, err = err)
    end function

! ----------
    subroutine sl_vecfcn(x, f, args)
        !! The subroutine passed to the inverse kinematics solver.
        real(real64), intent(in), dimension(:) :: x
            !! The N joint variables
        real(real64), intent(out), dimension(:) :: f
            !! The 6 kinematic constraint equations.
        class(*), intent(inout), optional :: args
            !! A container for the serial_linkage_solver_data object.

        ! Local Variables
        real(real64) :: T(4, 4), Re(3, 3), R(3, 3)

        ! Process
        select type (args)
        class is (serial_linkage_solver_data)
            ! Compute the forward kinematics of the linkage given the joint
            ! variables
            T = args%linkage%forward_kinematics(x)

            ! Evaluate the kinematics equations for position and orientation
            f(1:3) = T(1:3,4)

            ! The orientation components utilize an angle-axis approximation.
            ! The idea is to drive these 3 values to 0.
            Re = transpose(T(1:3,1:3))
            R = matmul(Re, args%target(1:3,1:3))
            f(4) = 0.5d0 * (R(3,2) - R(2,3))
            f(5) = 0.5d0 * (R(1,3) - R(3,1))
            f(6) = 0.5d0 * (R(2,1) - R(1,2))
        end select
    end subroutine

! ----------
    subroutine sl_jacobianfcn(x, jac, args)
        !! The Jacobian evaluation subroutine to pass to the solver.
        real(real64), intent(in), dimension(:) :: x
            !! The N joint variables
        real(real64), intent(out), dimension(:,:) :: jac
            !! The 6-by-N Jacobian.
        class(*), intent(inout), optional :: args
            !! A container for the serial_linkage_solver_data object.

        ! Process
        select type (args)
        class is (serial_linkage_solver_data)
            ! Compute the Jacobian matrix for the linkage
            jac = args%linkage%jacobian(x)
        end select
    end subroutine

! ------------------------------------------------------------------------------
end module