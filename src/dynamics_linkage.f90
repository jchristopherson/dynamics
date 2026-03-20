module dynamics_linkage
    use iso_fortran_env
    use dynamics_rigid_bodies
    use dynamics_kinematics
    use dynamics_geometry
    use dynamics_helper
    use dynamics_quaternions
    implicit none
    private
    public :: link
    public :: ground_link
    public :: initialize_ground_link
    public :: binary_link
    public :: initialize_binary_link

    type, abstract, extends(rigid_body) :: link
    contains
        procedure(link_integer_value), deferred, public, pass :: get_distal_joint_type
            !! Returns the joint type.  This value will be either REVOLUTE_JOINT
            !! or PRISMATIC_JOINT.
    end type

    interface
        pure function link_integer_value(this) result(rst)
            !! Defines the interface of a function for returning an integer
            !! parameter value from a link object.
            use iso_fortran_env, only : int32
            import link
            class(link), intent(in) :: this
                !! The link object.
            integer(int32) :: rst
                !! The result.
        end function
    end interface

    type, extends(link) :: ground_link
        !! Defines a ground link.
        integer(int32), private :: m_jointType
        real(real64), private :: m_jointAxis(3)
        real(real64), private :: m_linkXAxis(3)
    contains
        procedure, public :: get_distal_joint_type => gl_get_joint_type
        procedure, public :: get_joint_axis => gl_get_axis
        procedure, public :: get_link_x_axis => gl_get_x_axis
    end type

    interface ground_link
        module procedure :: gl_init
    end interface

    type, extends(link) :: binary_link
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
        integer(int32), private :: m_distalJoint
        integer(int32), private :: m_proximalJoint
    contains
        procedure, public :: get_proximal_joint_type => bl_get_proximal_joint
        procedure, public :: get_distal_joint_type => bl_get_distal_joint
    end type

    interface binary_link
        module procedure :: bl_init
    end interface

contains
! ******************************************************************************
! GROUND_LINK MEMBERS
! ------------------------------------------------------------------------------
    pure function gl_init(joint, axis, xaxis) result(rst)
        !! Initializes a new ground_link object.
        integer(int32), intent(in), optional :: joint
            !! The joint type.  This value must be either REVOLUTE_JOINT or
            !! PRISMATIC JOINT.  It defaults to REVOLUTE_JOINT.
        real(real64), intent(in), optional :: axis(3)
            !! The joint axis, as expressed in the parent coordinate frame.
            !! If not provided, this value defaults to [0, 0, 1].
        real(real64), intent(in), optional :: xaxis(3)
            !! The link x-axis as expressed in the parent coordinate frame.
            !! If not providied, this value defaults to [1, 0, 0].
        type(ground_link) :: rst
            !! The new ground_link object.

        call initialize_ground_link(rst, joint, axis, xaxis)
    end function

! ------------------------------------------------------------------------------
    pure subroutine initialize_ground_link(lnk, joint, axis, xaxis)
        !! Initializes a new ground_link object.
        class(ground_link), intent(inout) :: lnk
            !! The ground_link object.
        integer(int32), intent(in), optional :: joint
            !! The joint type.  This value must be either REVOLUTE_JOINT or
            !! PRISMATIC JOINT.  It defaults to REVOLUTE_JOINT.
        real(real64), intent(in), optional :: axis(3)
            !! The joint axis, as expressed in the parent coordinate frame.
            !! If not provided, this value defaults to [0, 0, 1].
        real(real64), intent(in), optional :: xaxis(3)
            !! The link x-axis as expressed in the parent coordinate frame.
            !! If not providied, this value defaults to [1, 0, 0].

        ! Initialization
        call initialize_rigid_body(lnk)
        if (present(joint)) then
            if (joint /= REVOLUTE_JOINT .and. joint /= PRISMATIC_JOINT) then
                lnk%m_jointType = REVOLUTE_JOINT
            else
                lnk%m_jointType = joint
            end if
        else
            lnk%m_jointType = REVOLUTE_JOINT
        end if
        if (present(axis)) then
            lnk%m_jointAxis = axis / norm2(axis)
        else
            lnk%m_jointAxis = [0.0d0, 0.0d0, 1.0d0]
        end if
        if (present(xaxis)) then
            lnk%m_linkXAxis = xaxis / norm2(xaxis)
        else
            lnk%m_linkXAxis = [1.0d0, 0.0d0, 0.0d0]
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure function gl_get_joint_type(this) result(rst)
        !! Gets the ground link joint type.  This value will be either 
        !! REVOLUTE_JOINT or PRISMATIC_JOINT.
        class(ground_link), intent(in) :: this
            !! The ground_link object.
        integer(int32) :: rst
            !! The joint type.
        rst = this%m_jointType
    end function

! ------------------------------------------------------------------------------
    pure function gl_get_axis(this) result(rst)
        !! Gets the axis of the joint in terms of the parent coordinate system.
        class(ground_link), intent(in) :: this
            !! The ground_link.
        real(real64) :: rst(3)
            !! A unit vector defining the joint axis.
        rst = this%m_jointAxis
    end function

! ------------------------------------------------------------------------------
    pure function gl_get_x_axis(this) result(rst)
        !! Gets the link x-axis in terms of the parent coordinate system.
        class(ground_link), intent(in) :: this
            !! The ground_link.
        real(real64) :: rst(3)
            !! A unit vector defining the link's x-axis.
        rst = this%m_linkXAxis
    end function

! ******************************************************************************
! BINARY_LINK MEMBERS
! ------------------------------------------------------------------------------
    pure function bl_init(proximal, distal, length, twist, offset, angle, &
        mass, inertia, cg) result(rst)
        !! Initializes a new parallel_revolute_revolute_link instance.
        integer(int32), intent(in), optional :: proximal
            !! The proximal joint type.  This value must be either &
            !! REVOLUTE_JOINT or PRISMATIC_JOINT.  If incorrectly specified, 
            !! this parameter defaults to REVOLUTE_JOINT.
        integer(int32), intent(in), optional :: distal
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

        ! Initialize the object
        call initialize_binary_link(rst, proximal, distal, length, twist, &
            offset, angle, mass, inertia, cg)
    end function

! ------------------------------------------------------------------------------
    pure subroutine initialize_binary_link(lnk, proximal, distal, length, &
        twist, offset, angle, mass, inertia, cg)
        !! Initializes a binary_link object.
        class(binary_link), intent(inout) :: lnk
            !! The binary_link object to initialize.
        integer(int32), intent(in), optional :: proximal
            !! The proximal joint type.  This value must be either &
            !! REVOLUTE_JOINT or PRISMATIC_JOINT.  If incorrectly specified, 
            !! this parameter defaults to REVOLUTE_JOINT.
        integer(int32), intent(in), optional :: distal
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

        ! Initialize the base object
        call initialize_rigid_body(lnk, mass, inertia, cg)
        
        ! Additional initialization
        if (present(proximal)) then
            if (proximal /= REVOLUTE_JOINT .and. proximal /= PRISMATIC_JOINT) then
                lnk%m_proximalJoint = REVOLUTE_JOINT
            else
                lnk%m_proximalJoint = proximal
            end if
        else
            lnk%m_proximalJoint = REVOLUTE_JOINT
        end if

        if (present(distal)) then
            if (distal /= REVOLUTE_JOINT .and. distal /= PRISMATIC_JOINT) then
                lnk%m_distalJoint = REVOLUTE_JOINT
            else
                lnk%m_distalJoint = distal
            end if
        else
            lnk%m_distalJoint = REVOLUTE_JOINT
        end if

        if (present(length)) then
            lnk%link_length = length
        else
            lnk%link_length = 0.0d0
        end if
        
        if (present(twist)) then
            lnk%link_twist = twist
        else
            lnk%link_twist = 0.0d0
        end if

        if (present(offset)) then
            lnk%link_offset = offset
        else
            lnk%link_offset = 0.0d0
        end if

        if (present(angle)) then
            lnk%joint_angle = angle
        else
            lnk%joint_angle = 0.0d0
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure function bl_get_proximal_joint(this) result(rst)
        !! Gets the type of the proximal joint.  This value will be either
        !! REVOLUTE_JOINT or PRISMATIC_JOINT.
        class(binary_link), intent(in) :: this
            !! The binary_link object.
        integer(int32) :: rst
            !! The joint type.
        rst = this%m_proximalJoint
    end function

! ------------------------------------------------------------------------------
    pure function bl_get_distal_joint(this) result(rst)
        !! Gets the type of the distal joint.  This value will be either
        !! REVOLUTE_JOINT or PRISMATIC_JOINT.
        class(binary_link), intent(in) :: this
            !! The binary_link object.
        integer(int32) :: rst
            !! The joint type.
        rst = this%m_distalJoint
    end function

! ------------------------------------------------------------------------------
end module