module dynamics_linkage
    use iso_fortran_env
    use dynamics_rigid_bodies
    use dynamics_kinematics
    implicit none

    type, extends(rigid_body) :: binary_link
        !! Defines a link consisting of only two joints.  The coordinate system
        !! of this link is situated at the distal joint with it's z-axis 
        !! coincident with the axis of the joint.
        integer(int32), public :: proximal_joint_type
            !! The proximal joint type.  This value must be either 
            !! PRISMATIC_JOINT or REVOLUTE_JOINT.
        integer(int32), public :: distal_joint_type
            !! The distal joint type.  This value must be either 
            !! PRISMATIC_JOINT or REVOLUTE_JOINT.
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
        real(real64), public :: angle
            !! The angle, in radians, required to bring the proximal joint's
            !! x-axis to parallel to the distal joint's x-axis as measured about
            !! the proximal joint's z-axis (proximal joint axis).
    end type

contains
! ******************************************************************************
! LINK MEMBERS & ROUTINES
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! KINEMATICS CALCULATIONS
! ------------------------------------------------------------------------------
pure function link_transformation_matrix(link, theta, d) result(rst)
    !! Computes the Denavit-Hartenberg transformation matrix for the specified
    !! binary_link object.  Assuming the link is labeled link \(i\), the 
    !! proximal joint is labeld \(i - 1 \) and is connects this link to the
    !! previous link (the \(i - 1\) link).  The computed transformation matrix
    !! is \(^{i-1}T_{i}) such that the following relationship holds.
    !! $$ \left[\begin{matrix} x_{i-1} \\ y_{i-1} \\ z_{i-1} \\ 1 = 
    !! \end{matrix}\right] = ^{i-1}T_{i} \left[\begin{matrix} x_{i} \\ y_{i} 
    !! \\ z_{i} \\ 1 \end{matrix}\right] $$
    class(binary_link), intent(in) :: link
        !! The binary_link object.
    real(real64), intent(in) :: theta
        !! The joint angle, in radians.
        !! specified.
    real(real64), intent(in) :: d
        !! The joint translation.
    real(real64) :: rst(4, 4)
        !! The resulting 4-by-4 link transformation matrix.

    ! Compute the matrix
    rst = dh_matrix(link%twist, link%length, link%angle + theta, link%offset + d)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module