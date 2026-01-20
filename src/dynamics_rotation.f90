module dynamics_rotation
    use iso_fortran_env
    use dynamics_helper
    implicit none
    private
    public :: rotate_x
    public :: rotate_y
    public :: rotate_z
    public :: rotate
    public :: acceleration_transform
    public :: velocity_transform
    public :: quaternion
    public :: abs
    public :: conjg
    public :: real
    public :: aimag

    interface rotate
        module procedure :: rotate_general_1
        module procedure :: rotate_general_2
    end interface

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
    end type

    interface quaternion
        module procedure :: quat_init_array
        module procedure :: quat_init_angle_axis
        module procedure :: quat_init_mtx
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
    
contains
! ------------------------------------------------------------------------------
pure function rotate_x(angle) result(rst)
        !! Constructs the rotation matrix describing a rotation about an
        !! x-axis such that 
        !! \( \overrightarrow{r_2} = \textbf{R}_x \overrightarrow{r_1} \).
        !!
        !! $$ \textbf{R}_x = \left[ \begin{matrix} 1 & 0 & 0 \\ 0 & 
        !! \cos{\theta_x} & -\sin{\theta_x} \\ 0 & \sin{\theta_x} & 
        !! \cos{\theta_x} \\ \end{matrix} \right] $$
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(3, 3)
            !! The resulting 3-by-3 matrix.

        ! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, c, s, 0.0d0, -s, c], [3, 3])
    end function

! ------------------------------------------------------------------------------
    pure function rotate_y(angle) result(rst)
        !! Constructs the rotation matrix describing a rotation about a y-axis
        !! such that 
        !! \( \overrightarrow{r_2} = \textbf{R}_y \overrightarrow{r_1} \).
        !!
        !! $$ \textbf{R}_y = \left[ \begin{matrix} \cos{\theta_y} & 0 & 
        !! \sin{\theta_y} \\ 0 & 1 & 0 \\ -\sin{\theta_y} & 0 & 
        !! \cos{\theta_y} \\ \end{matrix} \right] $$
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(3, 3)
            !! The resulting 3-by-3 matrix.

        ! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([c, 0.0d0, -s, 0.0d0, 1.0d0, 0.0d0, s, 0.0d0, c], [3, 3])
    end function

! ------------------------------------------------------------------------------
    pure function rotate_z(angle) result(rst)
        !! Constructs the rotation matrix describing a rotation about a y-axis
        !! such that 
        !! \( \overrightarrow{r_2} = \textbf{R}_z \overrightarrow{r_1} \).
        !!
        !! $$ \textbf{R}_z = \left[ \begin{matrix} \cos{\theta_z} & 
        !! -\sin{\theta_z} & 0 \\ \sin{\theta_z} & \cos{\theta_z} & 0 \\
        !! 0 & 0 & 1 \\ \end{matrix} \right] $$
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(3, 3)
            !! The resulting 3-by-3 matrix.

        ! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([c, s, 0.0d0, -s, c, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3, 3])
    end function

! ------------------------------------------------------------------------------
    pure function rotate_general_1(i, j, k, Ip, Jp, Kp) result(rst)
        !! Constructs a rotation matrix when the orientation of the coordinate
        !! frame of interest is known relative to the parent coordinate frame.
        !!
        !! The matrix is of the following form.
        !!
        !! $$ \textbf{R} = \left[ \begin{matrix} 
        !! \vec{I_p} \cdot \vec{i} & \vec{I_p} \cdot \vec{j} & 
        !! \vec{I_p} \cdot \vec{k} \\ \vec{J_p} \cdot \vec{i} &
        !! \vec{J_p} \cdot \vec{j} & \vec{J_p} \cdot \vec{k} \\ 
        !! \vec{K_p} \cdot \vec{i} & \vec{K_p} \cdot \vec{j} &
        !! \vec{K_p} \cdot \vec{k} \\ \end{matrix} \right] $$
        !!
        !! This routine does not check for orthogonallity or unit vector length;
        !! therefore, to ensure correct results it is the callers responsibility
        !! to ensure each vector is of unit length and that the unit vectors
        !! are properly orthogonal.
        real(real64), intent(in) :: i(3)
            !! The rotated coordinate frame x-axis unit vector.
        real(real64), intent(in) :: j(3)
            !! The rotated coordinate frame y-axis unit vector.
        real(real64), intent(in) :: k(3)
            !! The rotated coordinate frame z-axis unit vector.
        real(real64), intent(in) :: Ip(3)
            !! The parent coordinate frame x-axis unit vector.
        real(real64), intent(in) :: Jp(3)
            !! The parent coordinate frame y-axis unit vector.
        real(real64), intent(in) :: Kp(3)
            !! The parent coordinate frame z-axis unit vector.
        real(real64) :: rst(3, 3)
            !! The resulting 3-by-3 matrix.

        rst(1,1) = dot_product(Ip, i)
        rst(2,1) = dot_product(Jp, i)
        rst(3,1) = dot_product(Kp, i)

        rst(1,2) = dot_product(Ip, j)
        rst(2,2) = dot_product(Jp, j)
        rst(3,2) = dot_product(Kp, j)

        rst(1,3) = dot_product(Ip, k)
        rst(2,3) = dot_product(Jp, k)
        rst(3,3) = dot_product(Kp, k)
    end function

! ------------------------------------------------------------------------------
    pure function rotate_general_2(i, j, k) result(rst)
        !! Constructs a rotation matrix when the orientation of the coordinate
        !! frame of interest is known relative to the parent coordinate frame.
        !!
        !! The matrix is of the following form.
        !!
        !! $$ \textbf{R} = \left[ \begin{matrix} 
        !! \vec{I_p} \cdot \vec{i} & \vec{I_p} \cdot \vec{j} & 
        !! \vec{I_p} \cdot \vec{k} \\ \vec{J_p} \cdot \vec{i} &
        !! \vec{J_p} \cdot \vec{j} & \vec{J_p} \cdot \vec{k} \\ 
        !! \vec{K_p} \cdot \vec{i} & \vec{K_p} \cdot \vec{j} &
        !! \vec{K_p} \cdot \vec{k} \\ \end{matrix} \right] $$
        !!
        !! The parent coordinate frame is assumed to be as follows.
        !!
        !! $$ \vec{I_p} = \left( \begin{matrix} 1 & 0 & 0  \end{matrix} \right) $$
        !!
        !! $$ \vec{J_p} = \left( \begin{matrix} 0 & 1 & 0  \end{matrix} \right) $$
        !!
        !! $$ \vec{K_p} = \left( \begin{matrix} 0 & 0 & 1  \end{matrix} \right) $$
        !!
        !! This routine does not check for orthogonallity or unit vector length;
        !! therefore, to ensure correct results it is the callers responsibility
        !! to ensure each vector is of unit length and that the unit vectors
        !! are properly orthogonal.
        real(real64), intent(in) :: i(3)
            !! The rotated coordinate frame x-axis unit vector.
        real(real64), intent(in) :: j(3)
            !! The rotated coordinate frame y-axis unit vector.
        real(real64), intent(in) :: k(3)
            !! The rotated coordinate frame z-axis unit vector.
        real(real64) :: rst(3, 3)
            !! The resulting 3-by-3 matrix.

        rst = rotate_general_1(i, j, k, [1.0d0, 0.0d0, 0.0d0], &
            [0.0d0, 1.0d0, 0.0d0], [0.0d0, 0.0d0, 1.0d0])
    end function

! ******************************************************************************
! REVISION 1.0.8 ADDITIONS
! ------------------------------------------------------------------------------
    pure function acceleration_transform(alpha, omega, a, x) result(rst)
        !! Computes the acceleration transformation matrix relating the
        !! position of a point expressed in a rotating and translating body
        !! relative to its parent frame.
        !!
        !! The transformation matrix takes the following form.
        !!
        !! $$ A = \left[ \begin{matrix} \tilde{\alpha} - \tilde{\omega} 
        !! \tilde{\omega}^{T} & \vec{a} - \left( \tilde{\alpha} - \tilde{\omega} 
        !! \tilde{\omega}^{T} \right) \vec{x} \\ 0 & 0 \end{matrix} \right] $$
        !!
        !! where,
        !!
        !! $$ \tilde{\alpha} = \left[ \begin{matrix} 0 & -\alpha_z & \alpha_y \\
        !! \alpha_z & 0 & -\alpha_x \\ -\alpha_y & \alpha_x & 0 \end{matrix}
        !! \right] $$
        !!
        !! and,
        !!
        !! $$ \tilde{\omega} = \left[ \begin{matrix} 0 & -\omega_z & \omega_y \\
        !! \omega_z & 0 & -\omega_x \\ -\omega_y & \omega_x & 0 \end{matrix}
        !! \right] $$
        !!
        !! Given a vector describing the location on a moving body, 
        !! \(\vec{r_p}\), the matrix is used to report its acceleration 
        !! \(\vec{a_p} = A \vec{r_p}\).
        real(real64), intent(in) :: alpha(3)
            !! The angular acceleration vector.
        real(real64), intent(in) :: omega(3)
            !! The angular velocity vector.
        real(real64), intent(in) :: a(3)
            !! The translational acceleration vector describing the acceleration
            !! of the body in its parent coordinate frame.
        real(real64), intent(in) :: x(3)
            !! The position vector of the body in its parent coordinate frame.
        real(real64) :: rst(4, 4)
            !! The 4-by-4 transformation matrix.

        ! Compute alpha = alpha - omega * omega**T
        rst(1,1) = -omega(3)**2 - omega(2)**2
        rst(2,1) = omega(1) * omega(2) + alpha(3)
        rst(3,1) = omega(1) * omega(3) - alpha(2)
        rst(4,1) = 0.0d0

        rst(1,2) = omega(1) * omega(2) - alpha(3)
        rst(2,2) = -omega(3)**2 - omega(1)**2
        rst(3,2) = omega(2) * omega(3) + alpha(1)
        rst(4,2) = 0.0d0

        rst(1,3) = omega(1) * omega(3) + alpha(2)
        rst(2,3) = omega(2) * omega(3) - alpha(1)
        rst(3,3) = -omega(2)**2 - omega(1)**2
        rst(4,3) = 0.0d0
        
        ! Compute a - (alpha - omega * omega**T) x
        rst(1:3,4) = a - matmul(rst(1:3,1:3), x)
        rst(4,4) = 0.0d0
    end function

! ------------------------------------------------------------------------------
    pure function velocity_transform(omega, v, x) result(rst)
        !! Computes the velocity transformation matrix relating the position
        !! of a point expressed in a rotating and translating body relative to
        !! its parent frame.
        !!
        !! The transformation matrix takes the following form.
        !!
        !! $$ V = \left[ \begin{matrix} \tilde{\omega} & \vec{v} - 
        !! \tilde{\omega} \vec{x} \\ 0 & 0 \end{matrix} \right] $$
        !!
        !! where,
        !!
        !! $$ \tilde{\omega} = \left[ \begin{matrix} 0 & -\omega_z & \omega_y \\
        !! \omega_z & 0 & -\omega_x \\ -\omega_y & \omega_x & 0 \end{matrix}
        !! \right] $$
        !!
        !! Given a vector describing the location on a moving body, 
        !! \(\vec{r_p}\), the matrix is used to report its velocity 
        !! \(\vec{v_p} = V \vec{r_p}\).
        real(real64), intent(in) :: omega(3)
            !! The angular velocity vector.
        real(real64), intent(in) :: v(3)
            !! The translation velocity vector describing the velocity of the
            !! body in its parent coordinate frame.
        real(real64), intent(in) :: x(3)
            !! The position vector of the body in its parent coordinate frame.
        real(real64) :: rst(4, 4)
            !! The 4-by-4 transformation matrix.

        ! Process
        rst(1:3,1:3) = to_skew_symmetric(omega)
        rst(1:3,4) = v - matmul(rst(1:3,1:3), x)
        rst(4,:) = 0.0d0
    end function

! ******************************************************************************
! QUATERNIONS
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
        real(real64) :: esq(4), trc, e0, e1, e2, e3

        ! Process
        trc = r(1,1) + r(2,2) + r(3,3)
        esq(1) = 0.5d0 * (1.0d0 + trc)
        esq(2) = 0.25d0 * (1.0d0 + 2.0d0 * r(1,1) - trc)
        esq(3) = 0.25d0 * (1.0d0 + 2.0d0 * r(2,2) - trc)
        esq(4) = 0.2530 * (1.0d0 + 2.0d0 * r(3,3) - trc)

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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module