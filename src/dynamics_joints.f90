module dynamics_joints
    !! Provides types and routines describing the kinematic joints that connect
    !! the links of a mechanism.
    use iso_fortran_env
    use dynamics_error_handling
    use dynamics_kinematics, only : REVOLUTE_JOINT, PRISMATIC_JOINT
    use linalg, only : identity
    implicit none
    private
    public :: REVOLUTE_JOINT
    public :: PRISMATIC_JOINT
    public :: FIXED_JOINT
    public :: CYLINDRICAL_JOINT
    public :: UNIVERSAL_JOINT
    public :: SPHERICAL_JOINT
    public :: joint
    public :: joint_degrees_of_freedom

    integer(int32), parameter :: FIXED_JOINT = 2
        !! Defines a rigid connection permitting no relative motion.
    integer(int32), parameter :: CYLINDRICAL_JOINT = 3
        !! Defines a joint permitting rotation about, and translation along, the
        !! joint's z-axis.
    integer(int32), parameter :: UNIVERSAL_JOINT = 4
        !! Defines a joint permitting rotation about the joint's x-axis followed
        !! by rotation about the resulting y-axis.
    integer(int32), parameter :: SPHERICAL_JOINT = 5
        !! Defines a joint permitting rotation about all three axes.  The
        !! rotation is parameterized by a z-y-x Euler angle sequence.

    type joint
        !! Defines a kinematic joint connecting two links.  The joint attaches
        !! to a coordinate frame on each of the two links; the frames are
        !! identified by their index within each link.  The joint is coincident
        !! with both frames when all of its joint variables are zero.
        integer(int32), public :: joint_type = REVOLUTE_JOINT
            !! The type of joint.  This value must be one of FIXED_JOINT,
            !! REVOLUTE_JOINT, PRISMATIC_JOINT, CYLINDRICAL_JOINT,
            !! UNIVERSAL_JOINT, or SPHERICAL_JOINT.
        integer(int32), public :: parent_link = 0
            !! The index of the parent link.
        integer(int32), public :: parent_frame = 1
            !! The index of the coordinate frame on the parent link to which
            !! the joint attaches.
        integer(int32), public :: child_link = 0
            !! The index of the child link.
        integer(int32), public :: child_frame = 1
            !! The index of the coordinate frame on the child link to which the
            !! joint attaches.
        logical, public :: actuated = .false.
            !! True if the joint is driven; else, false for a passive joint.
    contains
        procedure, public :: get_dof => jnt_get_dof
        procedure, public :: get_constraint_count => jnt_get_constraint_count
        procedure, public :: get_planar_constraint_count => jnt_get_planar_count
        procedure, public :: is_planar_compatible => jnt_is_planar
        procedure, public :: motion_transform => jnt_motion_transform
    end type

    interface joint
        module procedure :: jnt_init
    end interface

contains
! ------------------------------------------------------------------------------
    function jnt_init(jtype, parent_link, child_link, parent_frame, &
        child_frame, actuated) result(rst)
        !! Initializes a new joint object.
        integer(int32), intent(in) :: jtype
            !! The type of joint.  This value must be one of FIXED_JOINT,
            !! REVOLUTE_JOINT, PRISMATIC_JOINT, CYLINDRICAL_JOINT,
            !! UNIVERSAL_JOINT, or SPHERICAL_JOINT.
        integer(int32), intent(in) :: parent_link
            !! The index of the parent link.
        integer(int32), intent(in) :: child_link
            !! The index of the child link.
        integer(int32), intent(in), optional :: parent_frame
            !! The index of the coordinate frame on the parent link to which the
            !! joint attaches.  If not supplied, a value of one is used.
        integer(int32), intent(in), optional :: child_frame
            !! The index of the coordinate frame on the child link to which the
            !! joint attaches.  If not supplied, a value of one is used.
        logical, intent(in), optional :: actuated
            !! True if the joint is driven.  If not supplied, the joint is
            !! considered passive.
        type(joint) :: rst
            !! The resulting joint object.

        ! Input Checking
        if (jtype /= FIXED_JOINT .and. jtype /= REVOLUTE_JOINT .and. &
            jtype /= PRISMATIC_JOINT .and. jtype /= CYLINDRICAL_JOINT .and. &
            jtype /= UNIVERSAL_JOINT .and. jtype /= SPHERICAL_JOINT) &
        then
            error stop DYN_INVALID_INPUT_ERROR
        end if
        if (parent_link < 1 .or. child_link < 1) then
            error stop DYN_INVALID_INPUT_ERROR
        end if

        ! Process
        rst%joint_type = jtype
        rst%parent_link = parent_link
        rst%child_link = child_link
        rst%parent_frame = 1
        rst%child_frame = 1
        rst%actuated = .false.
        if (present(parent_frame)) rst%parent_frame = parent_frame
        if (present(child_frame)) rst%child_frame = child_frame
        if (present(actuated)) rst%actuated = actuated
        if (rst%parent_frame < 1 .or. rst%child_frame < 1) then
            error stop DYN_INVALID_INPUT_ERROR
        end if
    end function

! ------------------------------------------------------------------------------
    pure function joint_degrees_of_freedom(jtype) result(rst)
        !! Gets the number of degrees of freedom permitted by the requested
        !! joint type.
        integer(int32), intent(in) :: jtype
            !! The joint type.
        integer(int32) :: rst
            !! The number of degrees of freedom.  A value of -1 is returned if
            !! the joint type is not recognized.

        select case (jtype)
        case (FIXED_JOINT)
            rst = 0
        case (REVOLUTE_JOINT, PRISMATIC_JOINT)
            rst = 1
        case (CYLINDRICAL_JOINT, UNIVERSAL_JOINT)
            rst = 2
        case (SPHERICAL_JOINT)
            rst = 3
        case default
            rst = -1
        end select
    end function

! ------------------------------------------------------------------------------
    pure function jnt_get_dof(this) result(rst)
        !! Gets the number of degrees of freedom, and therefore the number of
        !! joint variables, permitted by the joint.
        class(joint), intent(in) :: this
            !! The joint object.
        integer(int32) :: rst
            !! The number of degrees of freedom.

        rst = joint_degrees_of_freedom(this%joint_type)
    end function

! ------------------------------------------------------------------------------
    pure function jnt_get_constraint_count(this) result(rst)
        !! Gets the number of constraint equations the joint imposes upon the
        !! relative motion of the two links it connects when considered in
        !! three-dimensional space.
        class(joint), intent(in) :: this
            !! The joint object.
        integer(int32) :: rst
            !! The number of constraint equations.

        rst = 6 - this%get_dof()
    end function

! ------------------------------------------------------------------------------
    pure function jnt_get_planar_count(this) result(rst)
        !! Gets the number of constraint equations the joint imposes upon the
        !! relative motion of the two links it connects when the mechanism is
        !! restricted to planar motion.
        class(joint), intent(in) :: this
            !! The joint object.
        integer(int32) :: rst
            !! The number of constraint equations.

        rst = 3 - this%get_dof()
    end function

! ------------------------------------------------------------------------------
    pure function jnt_is_planar(this) result(rst)
        !! Determines if the joint is usable within a planar mechanism.  Planar
        !! mechanisms admit only fixed joints, revolute joints acting about the
        !! plane normal, and prismatic joints acting within the plane.
        class(joint), intent(in) :: this
            !! The joint object.
        logical :: rst
            !! True if the joint may be used in a planar mechanism.

        rst = this%joint_type == FIXED_JOINT .or. &
            this%joint_type == REVOLUTE_JOINT .or. &
            this%joint_type == PRISMATIC_JOINT
    end function

! ------------------------------------------------------------------------------
    function jnt_motion_transform(this, q) result(rst)
        !! Computes the transformation matrix describing the motion permitted by
        !! the joint.  The transformation is expressed in, and relative to, the
        !! joint's coordinate frame.
        !!
        !! For a planar mechanism the plane of motion is the x-y plane; a
        !! revolute joint therefore acts about the z-axis, and a prismatic joint
        !! acts along the joint frame's z-axis only when that axis lies in the
        !! plane of motion.
        class(joint), intent(in) :: this
            !! The joint object.
        real(real64), intent(in), dimension(:) :: q
            !! An array containing the joint variables.  The array must contain
            !! at least as many values as the joint has degrees of freedom.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Input Checking
        if (size(q) < this%get_dof()) error stop DYN_ARRAY_SIZE_ERROR

        ! Process
        select case (this%joint_type)
        case (FIXED_JOINT)
            rst = identity(4)
        case (REVOLUTE_JOINT)
            rst = rotate_z(q(1))
        case (PRISMATIC_JOINT)
            rst = translate_z(q(1))
        case (CYLINDRICAL_JOINT)
            rst = matmul(rotate_z(q(1)), translate_z(q(2)))
        case (UNIVERSAL_JOINT)
            rst = matmul(rotate_x(q(1)), rotate_y(q(2)))
        case (SPHERICAL_JOINT)
            rst = matmul(rotate_z(q(1)), matmul(rotate_y(q(2)), rotate_x(q(3))))
        case default
            error stop DYN_INVALID_INPUT_ERROR
        end select
    end function

! ------------------------------------------------------------------------------
    pure function rotate_x(angle) result(rst)
        ! A 4-by-4 rotation about the x-axis.
        real(real64), intent(in) :: angle
        real(real64) :: rst(4, 4)
        real(real64) :: c, s
        c = cos(angle)
        s = sin(angle)
        rst = identity(4)
        rst(2,2) = c
        rst(3,2) = s
        rst(2,3) = -s
        rst(3,3) = c
    end function

! ------------------------------------------------------------------------------
    pure function rotate_y(angle) result(rst)
        ! A 4-by-4 rotation about the y-axis.
        real(real64), intent(in) :: angle
        real(real64) :: rst(4, 4)
        real(real64) :: c, s
        c = cos(angle)
        s = sin(angle)
        rst = identity(4)
        rst(1,1) = c
        rst(3,1) = -s
        rst(1,3) = s
        rst(3,3) = c
    end function

! ------------------------------------------------------------------------------
    pure function rotate_z(angle) result(rst)
        ! A 4-by-4 rotation about the z-axis.
        real(real64), intent(in) :: angle
        real(real64) :: rst(4, 4)
        real(real64) :: c, s
        c = cos(angle)
        s = sin(angle)
        rst = identity(4)
        rst(1,1) = c
        rst(2,1) = s
        rst(1,2) = -s
        rst(2,2) = c
    end function

! ------------------------------------------------------------------------------
    pure function translate_z(d) result(rst)
        ! A 4-by-4 translation along the z-axis.
        real(real64), intent(in) :: d
        real(real64) :: rst(4, 4)
        rst = identity(4)
        rst(3,4) = d
    end function

! ------------------------------------------------------------------------------
end module
