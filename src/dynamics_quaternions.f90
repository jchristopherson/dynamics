module dynamics_quaternions
    use iso_fortran_env
    implicit none
    private
    public :: quaternion
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)
    public :: operator(/)
    public :: assignment(=)
    public :: operator(**)
    public :: abs
    public :: conjg
    public :: real
    public :: aimag
    public :: inverse
    public :: exp
    public :: log
    public :: dot_product

    type quaternion
        !! Defines a quaternion of the form: 
        real(real64) :: w
            !! The real component of the quaternion.
        real(real64) :: x
            !! The first element in the imaginary component of the quaternion. 
        real(real64) :: y
            !! The second element in the imaginary component of the quaternion.
        real(real64) :: z
            !! The third element in the imaginary component of the quaternion.
    contains
        procedure, public :: to_matrix => quat_to_matrix
        procedure, public :: normalize => quat_normalize
        procedure, public :: to_array => quat_to_vec4
        procedure, public :: to_angle_axis => quat_to_angle_axis
        procedure, public :: to_roll_pitch_yaw => quat_to_rpy
    end type

    interface quaternion
        module procedure :: quat_init_array
        module procedure :: quat_init_angle_axis
        module procedure :: quat_init_mtx
    end interface

    interface operator(+)
        module procedure :: quat_add
    end interface

    interface operator(-)
        module procedure :: quat_subtract
        module procedure :: quat_negate
    end interface

    interface operator(*)
        module procedure :: quat_scalar_mult
        module procedure :: scalar_quat_mult
        module procedure :: quat_multiply
        module procedure :: quat_vec3_mult
    end interface

    interface operator(/)
        module procedure :: quat_scalar_div
        module procedure :: quat_divide
    end interface

    interface assignment(=)
        module procedure :: quat_assign
        module procedure :: quat_assign_vec4
    end interface

    interface operator(**)
        module procedure :: quat_pwr
        module procedure :: quat_int_pwr
    end interface

    interface abs
        module procedure :: quat_abs
    end interface

    interface conjg
        module procedure :: quat_conjg
    end interface

    interface real
        module procedure :: quat_real
    end interface

    interface aimag
        module procedure :: quat_aimag
    end interface

    interface exp
        module procedure :: quat_exp
    end interface

    interface log
        module procedure :: quat_log
    end interface

    interface dot_product
        module procedure :: quat_dot_prd
    end interface

contains
! ------------------------------------------------------------------------------
    pure function quat_init_array(x) result(rst)
        !! Constructs a quaternion from a 4-element array stored such that
        !! [w, x, y, z].
        real(real64), intent(in), dimension(4) :: x
            !! The array from which to initialize the quaternion stored in the
            !! order [w, x, y, z].
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = x(1)
        rst%x = x(2)
        rst%y = x(3)
        rst%z = x(4)
    end function

! ------------------------------------------------------------------------------
    pure function quat_init_angle_axis(angle, axis) result(rst)
        !! Constructs a quaternion given an axis and the angle of rotation about
        !! the axis.
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64), intent(in), dimension(3) :: axis
            !! A 3-element vector defining the axis about which the rotation
            !! occurrs.
        type(quaternion) :: rst
            !! The resulting quaternion.

        real(real64) :: s
        
        s = sin(0.5d0 * angle)

        rst%w = cos(0.5d0 * angle)
        rst%x = s * axis(1)
        rst%y = s * axis(2)
        rst%z = s * axis(3)
    end function

! ------------------------------------------------------------------------------
    pure function quat_init_mtx(r) result(rst)
        !! Constructs a quaternion from a 3-by-3 rotation matrix using the 
        !! Stanley method.
        real(real64), intent(in), dimension(3, 3) :: r
            !! The rotation matrix.
        type(quaternion) :: rst
            !! The resulting quaternion.

        ! Local Variables
        integer(int32) :: ind
        real(real64) :: esq(4), e0, e1, e2, e3

        ! Process
        esq(1) = 0.25d0 * (1.0d0 + r(1,1) + r(2,2) + r(3,3))
        esq(2) = 0.25d0 * (1.0d0 + r(1,1) - r(2,2) - r(3,3))
        esq(3) = 0.25d0 * (1.0d0 - r(1,1) + r(2,2) - r(3,3))
        esq(4) = 0.25d0 * (1.0d0 - r(1,1) - r(2,2) + r(3,3))

        ind = maxloc(esq, 1)
        select case (ind)
        case (1)
            e0 = sqrt(esq(1))
            e1 = 0.25d0 * (r(3,2) - r(2,3)) / e0
            e2 = 0.25d0 * (r(1,3) - r(3,1)) / e0
            e3 = 0.25d0 * (r(2,1) - r(1,2)) / e0
        case (2)
            e1 = sqrt(esq(2))
            e0 = 0.25d0 * (r(3,2) - r(2,3)) / e1
            e2 = 0.25d0 * (r(1,2) + r(2,1)) / e1
            e3 = 0.25d0 * (r(1,3) + r(3,1)) / e1
        case (3)
            e2 = sqrt(esq(3))
            e0 = 0.25d0 * (r(1,3) - r(3,1)) / e2
            e1 = 0.25d0 * (r(1,2) + r(2,1)) / e2
            e3 = 0.25d0 * (r(2,3) + r(3,2)) / e2
        case default
            e3 = sqrt(esq(4))
            e0 = 0.25d0 * (r(2,1) - r(1,2)) / e3
            e1 = 0.25d0 * (r(1,3) + r(3,1)) / e3
            e2 = 0.25d0 * (r(2,3) + r(3,2)) / e3
        end select
        rst%w = e0
        rst%x = e1
        rst%y = e2
        rst%z = e3
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_add(x, y) result(rst)
        !! Adds two quaternions together.
        type(quaternion), intent(in) :: x
            !! The left-hand-side argument.
        type(quaternion), intent(in) :: y
            !! The right-hand-side argument.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = x%w + y%w
        rst%x = x%x + y%x
        rst%y = x%y + y%y
        rst%z = x%z + y%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_subtract(x, y) result(rst)
        !! Subtracts two quaternions.
        type(quaternion), intent(in) :: x
            !! The left-hand-side argument.
        type(quaternion), intent(in) :: y
            !! The right-hand-side argument.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = x%w - y%w
        rst%x = x%x - y%x
        rst%y = x%y - y%y
        rst%z = x%z - y%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_scalar_mult(x, y) result(rst)
        !! Multiplies a quaternion with a scalar.
        type(quaternion), intent(in) :: x
            !! The quaternion.
        real(real64), intent(in) :: y
            !! The scalar.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = y * x%w
        rst%x = y * x%x
        rst%y = y * x%y
        rst%z = y * x%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function scalar_quat_mult(x, y) result(rst)
        !! Multiplies a quaternion with a scalar.
        real(real64), intent(in) :: x
            !! The scalar.
        type(quaternion), intent(in) :: y
            !! The quaternion.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = x * y%w
        rst%x = x * y%x
        rst%y = x * y%y
        rst%z = x * y%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_multiply(x, y) result(rst)
        !! Multiplies two quaternions.
        type(quaternion), intent(in) :: x
            !! The left-hand-side argument.
        type(quaternion), intent(in) :: y
            !! The right-hand-side argument.
        type(quaternion) :: rst
            !! The resulting quaternion.
        
        ! Local Variables
        real(real64) :: x0, x1, x2, x3, y0, y1, y2, y3

        ! Initialization
        x0 = x%w
        x1 = x%x
        x2 = x%y
        x3 = x%z
        y0 = y%w
        y1 = y%x
        y2 = y%y
        y3 = y%z

        ! Process
        rst%w = x0 * y0 - x1 * y1 - x2 * y2 - x3 * y3
        rst%x = y0 * x1 + y1 * x0 - y2 * x3 + y3 * x2
        rst%y = y0 * x2 + x0 * y2 + y1 * x3 - x1 * y3
        rst%z = y0 * x3 - y1 * x2 + x0 * y3 + y2 * x1
    end function

! ------------------------------------------------------------------------------
    pure function quat_vec3_mult(x, y) result(rst)
        !! Multiplies a quaternion with a 3-element vector.
        type(quaternion), intent(in) :: x
            !! The quaternion.
        real(real64), intent(in), dimension(3) :: y
            !! The vector.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst = x * quaternion([0.0d0, y(1), y(2), y(3)])
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_scalar_div(x, y) result(rst)
        !! Divides a quaternion by a scalar.
        type(quaternion), intent(in) :: x
            !! The quaternion.
        real(real64), intent(in) :: y
            !! The scalar.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = x%w / y
        rst%x = x%x / y
        rst%y = x%y / y
        rst%z = x%z / y
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_divide(x, y) result(rst)
        !! Divides a quaternion by another.
        type(quaternion), intent(in) :: x
            !! The left-hand-side argument.
        type(quaternion), intent(in) :: y
            !! The right-hand-side argument.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst = x * inverse(y)
    end function

! ------------------------------------------------------------------------------
    pure elemental subroutine quat_assign(x, y)
        !! Assigns a quaternion to another.
        type(quaternion), intent(out) :: x
            !! The resulting quaternion.
        type(quaternion), intent(in) :: y
            !! The source quaternion.

        x%w = y%w
        x%x = y%x
        x%y = y%y
        x%z = y%z
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine quat_assign_vec4(x, y)
        !! Assigns a 4-element vector to a quaternion.
        type(quaternion), intent(out) :: x
            !! The resulting quaternion.
        real(real64), intent(in), dimension(4) :: y
            !! The source vector.

        x%w = y(1)
        x%x = y(2)
        x%y = y(3)
        x%z = y(4)
    end subroutine

! ------------------------------------------------------------------------------
    pure elemental function quat_negate(x) result(rst)
        !! Negates a quaternion.
        type(quaternion), intent(in) :: x
            !! The quaternion.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst%w = -x%w
        rst%x = -x%x
        rst%y = -x%y
        rst%z = -x%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_abs(q) result(rst)
        !! Computes the magnitude of a quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64) :: rst
            !! The magnitude of the quaternion.

        real(real64) :: x(4)

        x = [q%w, q%x, q%y, q%z]
        rst = norm2(x)
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_conjg(q) result(rst)
        !! Computes the conjugate of a quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        type(quaternion) :: rst
            !! The conjugate of the quaternion.

        rst%w = q%w
        rst%x = -q%x
        rst%y = -q%y
        rst%z = -q%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_real(q) result(rst)
        !! Returns the real component of the quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64) :: rst
            !! The real component.

        rst = q%w
    end function

! ------------------------------------------------------------------------------
    pure function quat_aimag(q) result(rst)
        !! Returns the imaginary component of the quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64), dimension(3) :: rst
            !! The imaginary component as a vector.

        rst = [q%x, q%y, q%z]
    end function

! ------------------------------------------------------------------------------
    pure elemental function inverse(q) result(rst)
        !! Computes the inverse of a quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        type(quaternion) :: rst
            !! The inverse of the supplied quaternion.

        rst = conjg(q) / (abs(q)**2)
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_exp(q) result(rst)
        !! Evaluates the exponential function for a quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        type(quaternion) :: rst
            !! The resulting quaternion.

        ! Local Variables
        real(real64) :: qmag, arg, s

        ! Process
        qmag = norm2(aimag(q))
        arg = exp(q%w)
        s = arg * sin(qmag) / qmag

        rst%w = arg * cos(qmag)
        rst%x = s * q%x
        rst%y = s * q%y
        rst%z = s * q%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_log(q) result(rst)
        !! Evalautes the natural logarithm function for a quaternion.
        type(quaternion), intent(in) :: q
            !! The quaternion.
        type(quaternion) :: rst
            !! The resulting quaternion.

        ! Local Variables
        real(real64) :: qmag, vmag, arg

        ! Process
        qmag = abs(q)
        vmag = norm2(aimag(q))
        arg = acos(q%w / qmag) / vmag

        rst%w = log(qmag)
        rst%x = arg * q%x
        rst%y = arg * q%y
        rst%z = arg * q%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_pwr(q, exponent) result(rst)
        !! Computes the result of a quaternion raised to an exponent.
        type(quaternion), intent(in) :: q
            !! The quaternion base.
        real(real64), intent(in) :: exponent
            !! The exponent.
        type(quaternion) :: rst
            !! The resulting quaternion.

        ! Local Variables
        real(real64) :: qmag, vmag, phi, arg, s

        ! Process
        qmag = abs(q)
        vmag = norm2(aimag(q))
        phi = acos(q%w / qmag)
        arg = qmag**exponent
        s = arg * sin(exponent * phi) / vmag

        rst%w = arg * cos(exponent * phi)
        rst%x = s * q%x
        rst%y = s * q%y
        rst%z = s * q%z
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_int_pwr(q, exponent) result(rst)
        !! Computes the result of a quaternion raised to an exponent.
        type(quaternion), intent(in) :: q
            !! The quaternion base.
        integer(int32), intent(in) :: exponent
            !! The exponent.
        type(quaternion) :: rst
            !! The resulting quaternion.

        rst = quat_pwr(q, real(exponent, real64))
    end function

! ------------------------------------------------------------------------------
    pure elemental function quat_dot_prd(x, y) result(rst)
        !! Computes the dot product of two quaternions.
        type(quaternion), intent(in) :: x
            !! The left-hand-side argument.
        type(quaternion), intent(in) :: y
            !! The right-hand-side argument.
        real(real64) :: rst
            !! The dot product.

        rst = x%w * y%w + x%x * y%x + x%y * y%y + x%z * y%z
    end function

! ******************************************************************************
! MEMBER ROUTINES
! ------------------------------------------------------------------------------
    pure function quat_to_matrix(q) result(rst)
        !! Converts the quaternion to a 3-by-3 rotation matrix.
        class(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64), dimension(3,3) :: rst
            !! The resulting rotation matrix.

        ! Local Variables
        type(quaternion) :: qnorm
        real(real64) :: e0, e1, e2, e3

        ! Process
        qnorm = q / abs(q)
        e0 = qnorm%w
        e1 = qnorm%x
        e2 = qnorm%y
        e3 = qnorm%z

        rst(1,1) = e0**2 + e1**2 - e2**2 - e3**2
        rst(2,1) = 2.0d0 * (e0 * e3 + e1 * e2)
        rst(3,1) = 2.0d0 * (e1 * e3 - e0 * e2)

        rst(1,2) = 2.0d0 * (e1 * e2 - e0 * e3)
        rst(2,2) = e0**2 - e1**2 + e2**2 - e3**2
        rst(3,2) = 2.0d0 * (e0 * e1 + e2 * e3)

        rst(1,3) = 2.0d0 * (e0 * e2 + e1 * e3)
        rst(2,3) = 2.0d0 * (e2 * e3 - e0 * e1)
        rst(3,3) = e0**2 - e1**2 - e2**2 + e3**2
    end function

! ------------------------------------------------------------------------------
    pure subroutine quat_normalize(q)
        !! Normalizes the quaternion.
        class(quaternion), intent(inout) :: q
            !! On input, the quaternion to normalize.  On output, the normalized
            !! quaternion.

        ! Local Variables
        real(real64) :: mag

        ! Process
        mag = abs(q)
        q%w = q%w / mag
        q%x = q%x / mag
        q%y = q%y / mag
        q%z = q%z / mag
    end subroutine

! ------------------------------------------------------------------------------
    pure function quat_to_vec4(q) result(rst)
        !! Converts a quaternion to a 4-element array of the form [w, x, y, z].
        class(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64), dimension(4) :: rst
            !! The array of the form [w, x, y, z].

        rst = [q%w, q%x, q%y, q%z]
    end function

! ------------------------------------------------------------------------------
    pure subroutine quat_to_angle_axis(q, angle, axis)
        !! Converts the quaternion to the equivalent angle-axis form.
        class(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64), intent(out) :: angle
            !! The rotation angle, in radians.
        real(real64), intent(out) :: axis(3)
            !! The axis of rotation.

        ! Process
        real(real64) :: qnorm, enorm, e0
        qnorm = abs(q)
        enorm = norm2(aimag(q))
        axis = aimag(q) / enorm
        angle = 2.0d0 * atan2(enorm, q%w / qnorm)
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine quat_to_rpy(q, roll, pitch, yaw)
        !! Converts the quaternion to the equivalent Euler angles of roll, 
        !! pitch, and yaw.
        class(quaternion), intent(in) :: q
            !! The quaternion.
        real(real64), intent(out) :: roll
            !! The roll angle (rotation about the body x-axis), in radians.
        real(real64), intent(out) :: pitch
            !! The pitch angle (rotation about the body y-axis), in radians.
        real(real64), intent(out) :: yaw
            !! The yaw angle (rotation about the body z-axis), in radians.

        ! Local Variables
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64) :: qmag, qx, qy, qz, qw

        ! Process
        qmag = abs(q)
        qw = q%w / qmag
        qx = q%x / qmag
        qy = q%y / qmag
        qz = q%z / qmag

        roll = atan2(2.0d0 * (qw * qx + qy * qz), 1.0d0 - 2.0d0 * (qx**2 + qy**2))
        pitch = -0.5d0 * pi + 2.0d0 * &
            atan2(sqrt(1.0d0 + 2.0d0 * (qw * qy - qx * qz)), &
            sqrt(1.0d0 - 2.0d0 * (qw * qy - qx * qz)))
        yaw = atan2(2.0d0 * (qw * qz + qx * qy), 1.0d0 - 2.0d0 * (qy**2 + qz**2))
    end subroutine

! ------------------------------------------------------------------------------
end module