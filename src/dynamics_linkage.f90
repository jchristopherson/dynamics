module dynamics_linkage
    use iso_fortran_env
    use dynamics_rigid_bodies
    use dynamics_kinematics
    use dynamics_geometry
    use dynamics_helper
    use dynamics_quaternions
    implicit none
    private

    type, abstract, extends(rigid_body) :: link
    end type


    type, extends(rigid_body) :: dh_binary_link
        !! Defines a link consisting of only two joints.  The coordinate system
        !! of this link is situated at the distal joint with it's z-axis 
        !! coincident with the axis of the joint.  The link utilizes a
        !! Denavit-Hartenberg convention in order to express its geometry.
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

    interface dh_binary_link
        module procedure :: dbl_init_1
    end interface

contains

! ******************************************************************************
! LINK MEMBERS & ROUTINES
! ------------------------------------------------------------------------------
pure function dbl_init_1(zim1, zi,  xim1, rim1, ri, jim1, ji) result(rst)
    !! Constructs a new dh_binary_link object.
    real(real64), intent(in) :: zim1(3)
        !! A 3-element vector describing the orientation of the proximal joint
        !! axis in terms of an external (parent) coordinate system.
    real(real64), intent(in) :: zi(3)
        !! A 3-element vector describing the orientation of the distal joint
        !! axis in terms of an external (parent) coordinate system.
    real(real64), intent(in) :: xim1(3)
        !! A 3-element vector describing the orientation of the prior link's 
        !! x-axis in terms of an external (parent) coordinate system.
    real(real64), intent(in) :: rim1(3)
        !! A 3-element vector describing the position of the proximal joint
        !! in terms of an external (parent) coordinate system.
    real(real64), intent(in) :: ri(3)
        !! A 3-element vector describing the position of the distal joint
        !! in terms of an external (parent) coordinate system.
    integer(int32), intent(in), optional :: jim1
        !! The type of joint for the proximal joint.  This value must be either
        !! PRISMATIC_JOINT or REVOLUTE_JOINT.  The joint type defaults to
        !! REVOLUTE_JOINT if not specified or specified incorrectly.
    integer(int32), intent(in), optional :: ji
        !! The type of joint for the distal joint.  This value must be either
        !! PRISMATIC_JOINT or REVOLUTE_JOINT.  The joint type defaults to
        !! REVOLUTE_JOINT if not specified or specified incorrectly.
    type(dh_binary_link) :: rst
        !! The resulting dh_binary_link object.

    ! Local Variables
    type(line) :: lzi, lzim1, lxi, lyi
    real(real64) :: tol, pzim1(3)

    ! Initialization
    tol = 1.0d1 * epsilon(tol)
    lzi = line_from_point_and_vector(ri, zi, .true.)
    lzim1 = line_from_point_and_vector(rim1, zim1, .true.)

    ! -------------------- LINK LENGTH -------------------- !
    ! Define the link's x-axis (points from proximal to distal)
    lxi = line_common_normal(lzim1, lzi)
    rst%length = norm2(lxi%evaluate(1.0d0) - lxi%evaluate(0.0d0)) ! link length
    lxi%r0 = ri ! ensure the x_i axis is correctly located

    ! If the two axes intersect, line_common_normal returns a zero-length line
    ! We then need to define the actual length
    if (rst%length < tol) then
        lxi%v = cross_product(lzi%v, lzim1%v)
    end if
    lxi%v = lxi%v / norm2(lxi%v)

    ! -------------------- LINK TWIST -------------------- !

    ! Project the z_i-1 axis onto the plane formed by the z_i and y_i axes.
    ! From this configuration, we can easily determine the rotation angle 
    ! required to bring the z_i-1 axis to parallel with the z_i axis.  The
    ! result is the link twist angle (alpha)
    lyi = line_from_point_and_vector(ri, cross_product(lzi%v, lxi%v), .true.)

    if (is_parallel(lxi%v, zim1)) then
        rst%twist = 0.0d0
    else
        ! To project a vector onto a plane, assume we're projecting vector X 
        ! onto plane P that has a unit normal N.  The projected vector is 
        ! X - (X dot N) * N.  For this instance, the normal vector of the y-z 
        ! plane is the x_i axis.  We've already normalized this vector to a 
        ! unit vector; therefore, the projected vector is computed as follows.
        pzim1 = zim1 - dot_product(zim1, lxi%v) * lxi%v

        ! Now we can compute the twist about the x_i axis by computing the angle
        ! between z_i and the projected z_i-1 axes
        rst%twist = vector_angle(lzi%v, pzim1)

        ! TO DO: Determine the proper sign for rst%twist
    end if

    ! -------------------- LINK OFFSET -------------------- !

    ! -------------------- JOINT ANGLE -------------------- !

    ! -------------------- INERTIAL PROPERTIES -------------------- !
end function

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
    !! is such that the following relationship holds.
    !!
    !! $$ \left[\begin{matrix} x_{i-1} \\ y_{i-1} \\ z_{i-1} \\ 1 
    !! \end{matrix}\right] = ^{i-1}T_{i} \left[\begin{matrix} x_{i} \\ y_{i} 
    !! \\ z_{i} \\ 1 \end{matrix}\right] $$
    !!
    !! In another words, this transformation matrix allows expression of a
    !! vector expressed in this link's coordinate system in terms of the 
    !! previous link's coordinate system.
    class(dh_binary_link), intent(in) :: link
        !! The dh_binary_link object.
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