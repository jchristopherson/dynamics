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


    type, abstract, extends(link) :: binary_link
        !! Defines a link consisting of only two joints.  The coordinate system
        !! of this link is situated at the distal joint with it's z-axis 
        !! coincident with the axis of the joint.  The link utilizes a
        !! Denavit-Hartenberg convention in order to express its geometry.   
        real(real64), public :: length
            !! The link length measured as the distance between the proximal and
            !! distal joints as projected along the link's x-axis.
        real(real64), public :: twist
            !! The twist angle of the link, in radians, that is required to 
            !! bring the proximal joint axis to parallel with the distal joint
            !! axis as measured about the link's x-axis.
        real(real64), public :: offset
            !! The distance between the proximal and distal joint's x-axes as
            !! measured along the proximal joint's z-axis.
    contains
        procedure(bl_get_integer), deferred, public, pass :: get_proximal_joint_type
            !! The proximal joint type.  This value will be either 
            !! PRISMATIC_JOINT or REVOLUTE_JOINT.
        procedure(bl_get_integer), deferred, public, pass :: get_distal_joint_type
            !! The distal joint type.  This value must be either 
            !! PRISMATIC_JOINT or REVOLUTE_JOINT.
    end type

    interface
        pure function bl_get_integer(this) result(rst)
            !! Defines the interface for a routine for retrieving integer 
            !! property values from a binary_link object.
            use iso_fortran_env, only : int32
            import binary_link
            class(binary_link), intent(in) :: this
                !! The binary_link object.
            integer(int32) :: rst
                !! The integer value.
        end function
    end interface

contains

! ******************************************************************************
! LINK MEMBERS & ROUTINES
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module