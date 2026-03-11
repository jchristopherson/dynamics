module dynamics_linkage
    use iso_fortran_env
    use dynamics_rigid_bodies
    use dynamics_kinematics
    use dynamics_geometry
    use dynamics_helper
    use dynamics_quaternions
    implicit none
    private
    public :: coordinate_system
    public :: dh_parameter_set
    public :: dh_table

    type coordinate_system
        !! Defines a 3D Cartesian coordinate system.
        real(real64) :: origin(3)
            !! The location of the origin of the coordinate system in the parent
            !! coordinate system.
        real(real64) :: i(3)
            !! The unit vector along the coordinate system's x-axis.
        real(real64) :: j(3)
            !! The unit vector along the coordinate system's y-axis.
        real(real64) :: k(3)
            !! The unit vector along the coordinate system's z-axis.
    end type

    interface coordinate_system
        module procedure :: define_link_csys
        module procedure :: define_csys
    end interface

    type dh_parameter_set
        !! Describes a set of Denavit-Hartenberg parameters for a single joint.
        real(real64) :: link_length
            !! The link length is the distance between the proximal and distal
            !! joint axes as measured along the link's x-axis.
        real(real64) :: link_twist
            !! The link twist is the required rotation of the proximal joint 
            !! axis about the link's x-axis to become parallel to the distal
            !! joint's axis.
        real(real64) :: link_offset
            !! The link offset is the distance between the previous link's
            !! x-axis and the current link's x-axis as measured along the
            !! axis of the proximal joint.
        real(real64) :: joint_angle
            !! The joint angle is the required rotation of the previous link's
            !! x-axis about the proximal joint's axis to become parallel to the
            !! current link's x-axis.
    end type

    interface dh_parameter_set
        module procedure :: dps_init
    end interface

    type dh_table
        !! Describes a Denavit-Hartenberg table.
        type(dh_parameter_set), allocatable, dimension(:) :: parameters
            !! A collection of DH parameters for each joint in the linkage.
    end type

    interface dh_table
        module procedure :: dh_table_init
    end interface

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
! COORDINATE SYSTEM ROUTINES
! ------------------------------------------------------------------------------
    pure function define_link_csys(xim1, zim1, zi, rim1, ri) &
        result(rst)
        !! Defines the DH coordinate system for the specified link.
        real(real64), intent(in) :: xim1(3)
            !! The x-axis of the previous link.
        real(real64), intent(in) :: zim1(3)
            !! The z-axis of the previous link.  This is also the axis of the
            !! proximal joint of the current link.
        real(real64), intent(in) :: zi(3)
            !! The axis of the distal joint of the current link.
        real(real64), intent(in) :: rim1(3)
            !! The location of the proximal joint's center or home position.
        real(real64), intent(in) :: ri(3)
            !! The location of the distal joint's center or home position.
        type(coordinate_system) :: rst
            !! The resulting coordinate system.

        ! Local Variables
        logical :: intersect
        type(line) :: cn, lzim1, lzi
        real(real64) :: length, tol, pt0(3)

        ! Initialization
        tol = 1.0d1 * epsilon(tol)  ! zero tolerance
        lzim1 = line_from_point_and_vector(rim1, zim1)
        lzi = line_from_point_and_vector(ri, zi)

        ! Process
        call do_lines_intersect(lzim1, lzi, intersect)
        if (intersect) then
            ! The two joint axes intersect.  The x-axis as follows.
            rst%i = cross_product(zim1, zi)

            ! Define the origin
            rst%origin = ri
        else
            ! Define the common normal between the two axes
            cn = line_common_normal(lzim1, lzi)
            pt0 = cn%evaluate(1.0d0)
            length = norm2(pt0 - cn%evaluate(0.0d0))
            if (length < tol) then
                ! The axes are effectively colinear.  The x-axis is chosen
                ! such that theta = 0 (i.e. the x axis is parallel with the
                ! previous x axis)
                rst%i = xim1
                pt0 = ri
            else
                rst%i = cn%v
            end if
            
            ! The origin is the point of intersection of the common normal with
            ! the distal joint axis
            rst%origin = pt0
        end if
        rst%i = rst%i / norm2(rst%i)    ! ensure i is a unit vector

        ! Store z, and normalize to a unit vector
        rst%k = zi / norm2(zi)

        ! Compute y
        rst%j = cross_product(rst%k, rst%i)
    end function

! ------------------------------------------------------------------------------
    pure function define_csys(i, j, k, o) result(rst)
        !! Defines a new coordinate_system object.  It is the callers 
        !! responsibility to ensure that the supplied vectors are orthogonal
        !! to one another.
        real(real64), intent(in) :: i(3)
            !! The x-axis unit vector.
        real(real64), intent(in) :: j(3)
            !! The y-axis unit vector.
        real(real64), intent(in) :: k(3)
            !! The x-axis unit vector.
        real(real64), intent(in), optional :: o(3)
            !! The location of the coordinate system within a reference 
            !! coordinate system.  If not supplied, a value of (0, 0, 0) will
            !! be used.
        type(coordinate_system) :: rst
            !! The new coordinate_system object.

        rst%i = i / norm2(i)
        rst%j = j / norm2(j)
        rst%k = k / norm2(k)
        if (present(o)) then
            rst%origin = o
        else
            rst%origin = [0.0d0, 0.0d0, 0.0d0]
        end if
    end function

! ******************************************************************************
! DENAVIT-HARTENBERG PARAMETER SET ROUTINES
! ------------------------------------------------------------------------------
    pure function dps_init(length, twist, offset, angle) result(rst)
        !! Constructs a new dh_parameter_set object.
        real(real64), intent(in) :: length
            !! The link length.
        real(real64), intent(in) :: twist
            !! The link twist.
        real(real64), intent(in) :: offset
            !! The link offset.
        real(real64), intent(in) :: angle
            !! The joint angle.
        type(dh_parameter_set) :: rst
            !! The new dh_parameter_set object.

        rst%link_length = length
        rst%link_twist = twist
        rst%link_offset = offset
        rst%joint_angle = angle
    end function

! ******************************************************************************
! DENAVIT-HARTENBERG TABLE ROUTINES
! ------------------------------------------------------------------------------
    pure function dh_table_init(c) result(rst)
        !! Constructs a Denavit-Hartenberg table given the locations and 
        !! orientations of a linkage.
        class(coordinate_system), intent(in), dimension(:) :: c
            !! An N-element array containing the coordinate frames for each
            !! of the N links of the linkage.
        type(dh_table) :: rst
            !! The resulting Denavit-Hartenberg table.

        ! Local Variables
        integer(int32) :: i, n

        ! Initialization
        n = size(c)
        allocate(rst%parameters(n - 1))

        ! Process
        do i = 2, n
            rst%parameters(i - 1)%link_length = &
                compute_dh_link_length(c(i - 1), c(i))
            rst%parameters(i - 1)%link_twist = &
                compute_dh_link_twist(c(i - 1), c(i))
            rst%parameters(i - 1)%link_offset = &
                compute_dh_link_offset(c(i - 1), c(i))
            rst%parameters(i - 1)%joint_angle = &
                compute_dh_joint_angle(c(i - 1), c(i))
        end do
    end function

! ------------------------------------------------------------------------------
    pure function compute_dh_link_length(cs1, cs2) result(rst)
        !! Computes the Denavit-Hartenberg link length as the distance between
        !! the cs2 z-axis and the cs1 z-axis as measured along the cs2 x-axis.
        class(coordinate_system), intent(in) :: cs1
            !! The previous coordinate frame.
        class(coordinate_system), intent(in) :: cs2
            !! The current coordinate frame.
        real(real64) :: rst
            !! The link length.

        ! Local Variables
        real(real64) :: v(3), pv(3)

        ! Compute a vector between the two coordinate system origins
        v = cs2%origin - cs1%origin

        ! Project the vector onto the cs2 x-axis
        pv = vector_projection(v, cs2%i)

        ! The result is the length of the projected vector
        rst = norm2(pv)
    end function

! ------------------------------------------------------------------------------
    pure function compute_dh_link_twist(cs1, cs2) result(rst)
        !! Computes the Denavit-Hartenberg link twist as the required rotation
        !! of the cs1 z-axis about the cs2 x-axis to become parallel to the
        !! cs2 z-axis.
        class(coordinate_system), intent(in) :: cs1
            !! The previous coordinate frame.
        class(coordinate_system), intent(in) :: cs2
            !! The current coordinate frame.
        real(real64) :: rst
            !! The link twist angle, in radians.  A positive angle is defined
            !! via the right-hand-rule.

        ! Local Variables
        real(real64) :: pz1(3), tol

        ! Initialization
        tol = 1.0d1 * epsilon(tol)

        ! Description:
        ! To compute the angle about the cs2 x-axis, the cs1 z-axis will be
        ! projected onto the plane formed by the cs2 y and z axes.  The angle
        ! between the projected vector and the cs2 z-axis can then be 
        ! determined.
        !
        ! To project a vector onto a plane, assume we're projecting vector X 
        ! onto plane P that has a unit normal N.  The projected vector is then
        ! X - (X dot N) * N.  For this instance, the normal vector of the cs2
        ! y-z plane is the cs2 x-axis.
        pz1 = cs1%k - dot_product(cs1%k, cs2%i) * cs2%i

        ! Now, the angle between vectors may be determined.
        if (norm2(pz1) < tol) then
            ! The projected vector has a length of zero
            rst = 0.0d0
        else
            ! Compute the angle
            rst = compute_vector_angle(cs2%i, cs2%k, pz1)
        end if
    end function

! ------------------------------------------------------------------------------
    pure function compute_dh_link_offset(cs1, cs2) result(rst)
        !! Computes the Denavit-Hartenberg link offset as the distance between
        !! the cs1 x-axis and the cs2 x-axis as measured along the cs1 z-axis.
        class(coordinate_system), intent(in) :: cs1
            !! The previous coordinate frame.
        class(coordinate_system), intent(in) :: cs2
            !! The current coordinate frame.
        real(real64) :: rst
            !! The link offset.

        ! Local Variables
        real(real64) :: v(3), pv(3)

        ! Compute a vector between the two coordinate system origins
        v = cs2%origin - cs1%origin

        ! Project the vector onto the cs1 z-axis
        pv = vector_projection(v, cs1%k)

        ! The result is the length of the projected vector
        rst = norm2(pv)
    end function

! ------------------------------------------------------------------------------
    pure function compute_dh_joint_angle(cs1, cs2) result(rst)
        !! Computes the Denavit-Hartenberg joint angle as rotation required of 
        !! cs1 x-axis about the cs1 z-axis to become parallel to the cs2 x-axis.
        class(coordinate_system), intent(in) :: cs1
            !! The previous coordinate frame.
        class(coordinate_system), intent(in) :: cs2
            !! The current coordinate frame.
        real(real64) :: rst
            !! The joint angle, in radians.  A positive angle is defined
            !! via the right-hand-rule.

        ! Local Variables
        real(real64) :: px2(3), tol

        ! Initialization
        tol = 1.0d1 * epsilon(tol)

        ! Description:
        ! To compute the angle about the cs1 z-axis, the cs2 x-axis will be
        ! projected onto the plane formed by the cs1 x-y plane.  The angle 
        ! between the projected vector and the cs1 x-axis can then be 
        ! determined.
        !
        ! To project a vector onto a plane, assume we're projecting vector X 
        ! onto plane P that has a unit normal N.  The projected vector is then
        ! X - (X dot N) * N.  For this instance, the normal vector of the cs1
        ! x-y plane is the cs1 z-axis.
        px2 = cs2%i - dot_product(cs2%i, cs1%k) * cs1%k

        ! Now, the angle between vectors may be determined
        if (norm2(px2) < tol) then
            ! The projected vector has a length of zero
            rst = 0.0d0
        else
            ! Compute the angle
            rst = compute_vector_angle(cs1%k, cs1%i, px2)
        end if
    end function

! ------------------------------------------------------------------------------
    pure function compute_vector_angle(axis, x, y) result(rst)
        !! Computes the angle of rotation between two vectors as measured about
        !! the specified axis.
        real(real64), intent(in) :: axis(3)
            !! The rotation axis.
        real(real64), intent(in) :: x(3)
            !! The target vector.
        real(real64), intent(in) :: y(3)
            !! The other vector.
        real(real64) :: rst
            !! The angle, with correct sign assuming a right-hand convention.

        ! Local Variables
        type(quaternion) :: q
        real(real64) :: ry(3)

        ! Determine the angle
        rst = vector_angle(x, y)

        ! To determine the sign, we can establish a quaternion and perform the
        ! rotation of one vector onto the other.  If the two vectors become
        ! parallel, we guessed the sign correctly.  If not, the sign is the
        ! opposite.
        q = quaternion(-rst, axis)   ! angle and axis construction
        ry = aimag(q * y * conjg(q))
        if (.not.is_parallel(ry, x)) rst = -rst
    end function

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