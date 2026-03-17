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
    public :: binary_link

    type, abstract, extends(rigid_body) :: link
    end type

    type, extends(link) :: binary_link
        !! Defines a link consisting of only two joints.  The coordinate system
        !! of this link is situated at the distal joint with it's z-axis 
        !! coincident with the axis of the joint.  The link utilizes a
        !! Denavit-Hartenberg convention in order to express its geometry.
        real(real64) :: link_length
            !! The link length is the distance between the proximal and distal
            !! joint axes as measured along the link's x-axis.
        real(real64) :: link_twist
            !! The link twist is the required rotation of the proximal joint 
            !! axis about the link's x-axis to become parallel to the distal
            !! joint's axis.
        real(real64) :: link_offset
            !! The link offset is the fixed distance between the previous link's
            !! x-axis and the current link's x-axis as measured along the
            !! axis of the proximal joint.
        real(real64) :: joint_angle
            !! The joint angle is the required rotation of the previous link's
            !! x-axis about the proximal joint's axis to become parallel to the
            !! current link's x-axis.
    end type

    interface binary_link
        module procedure :: bl_init
    end interface

contains
! ******************************************************************************
! BINARY_LINK MEMBERS
! ------------------------------------------------------------------------------
    pure function bl_init(length, twist, offset, angle, mass, inertia, cg) result(rst)
        !! Initializes a new parallel_revolute_revolute_link instance.
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

! ------------------------------------------------------------------------------
end module