! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
module dynamics_rotation
    use iso_fortran_env
    use dynamics_helper
    use linalg, only : identity
    implicit none
    private
    public :: rotate_x
    public :: rotate_y
    public :: rotate_z
    public :: homogeneous_rotation_x
    public :: homogeneous_rotation_y
    public :: homogeneous_rotation_z
    public :: rotate
    public :: translate
    public :: acceleration_transform
    public :: velocity_transform
    public :: to_angle_axis

    interface rotate
        module procedure :: rotate_general_1
        module procedure :: rotate_general_2
    end interface

    interface translate
        module procedure :: translate_1
        module procedure :: translate_2
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
    pure function homogeneous_rotation_x(angle) result(rst)
        !! Constructs the 4-by-4 homogeneous transformation matrix describing
        !! a rotation about an x-axis.
        !!
        !! $$ \textbf{R}_x = \left[ \begin{matrix} 1 & 0 & 0 & 0 \\ 0 &
        !! \cos{\theta_x} & -\sin{\theta_x} & 0 \\ 0 & \sin{\theta_x} &
        !! \cos{\theta_x} & 0 \\ 0 & 0 & 0 & 1 \\ \right] $$
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        !! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([ &
            1.0d0, 0.0d0, 0.0d0, 0.0d0, &
            0.0d0, c, s, 0.0d0, &
            0.0d0, -s, c, 0.0d0, &
            0.0d0, 0.0d0, 0.0d0, 1.0d0], &
            [4, 4])
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
    pure function homogeneous_rotation_y(angle) result(rst)
        !! Constructs the 4-by-4 homogeneous transformation matrix describing
        !! a rotation about a y-axis.
        !!
        !! $$ \textbf{R}_y = \left[ \begin{matrix} \cos{\theta_y} & 0 & 
        !! \sin{\theta_y} & 0 \\ 0 & 1 & 0 & 0 \\ -\sin{\theta_y} & 0 & 
        !! \cos{\theta_y} & 0 \\ 0 & 0 & 0 & 1 \\ \end{matrix} \right] $$
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 matrix.

        ! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([ &
            c, 0.0d0, -s, 0.0d0, &
            0.0d0, 1.0d0, 0.0d0, 0.0d0, &
            s, 0.0d0, c, 0.0d0, &
            0.0d0, 0.0d0, 0.0d0, 1.0d0], &
            [4, 4])
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
    pure function homogeneous_rotation_z(angle) result(rst)
        !! Constructs the 4-by-4 homogeneous transformation matrix describing
        !! a rotation about a y-axis.
        !!
        !! $$ \textbf{R}_z = \left[ \begin{matrix} \cos{\theta_z} & 
        !! -\sin{\theta_z} & 0 & 0 \\ \sin{\theta_z} & \cos{\theta_z} & 0 & 0 \\
        !! 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ \end{matrix} \right] $$
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 matrix.

        ! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([c, s, 0.0d0, 0.0d0, &
            -s, c, 0.0d0, 0.0d0, &
            0.0d0, 0.0d0, 1.0d0, 0.0d0, &
            0.0d0, 0.0d0, 0.0d0, 1.0d0], &
            [4, 4])
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
        !! For orthonormal frames this is a proper rotation, so
        !! \(R^T R=I\) and \(\det(R)=1\).
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

! ------------------------------------------------------------------------------
    pure function translate_1(x, y, z) result(rst)
        !! Computes the 4-by-4 homogeneous transformation matrix describing a
        !! rigid-body translation.
        real(real64), intent(in) :: x
            !! The x-component of the translation.
        real(real64), intent(in) :: y
            !! The y-component of the translation.
        real(real64), intent(in) :: z
            !! The z-component of the translation.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 matrix.

        rst = identity(4)
        rst(1,4) = x
        rst(2,4) = y
        rst(3,4) = z
    end function

! ------------------------------------------------------------------------------
    pure function translate_2(d) result(rst)
        !! Computes the 4-by-4 homogeneous transformation matrix describing a
        !! rigid-body translation.
        real(real64), intent(in) :: d(3)
            !! The x, y, z translation vector.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 matrix.

        rst = identity(4)
        rst(1,4) = d(1)
        rst(2,4) = d(2)
        rst(3,4) = d(3)
    end function

! ------------------------------------------------------------------------------
    pure subroutine to_angle_axis(r, angle, axis)
        !! Extracts the equivalent rotation angle and axis of rotation given a
        !! 3-by-3 rotation matrix.
        !! For a proper rotation, the angle is recovered from
        !! $$ \theta = \cos^{-1}\left(\frac{\operatorname{tr}(R)-1}{2}\right), $$
        !! while the skew-symmetric part satisfies
        !! $$ R-R^T = 2\sin(\theta)[\hat{u}]_\times. $$
        real(real64), intent(in) :: r(3,3)
            !! The 3-by-3 rotation matrix.
        real(real64), intent(out) :: angle
            !! The rotation angle.
        real(real64), intent(out) :: axis(3)
            !! The axis of rotation.

        ! Process
        real(real64) :: u(3,3)
        angle = acos(0.5d0 * (r(1,1) + r(2,2) + r(3,3) - 1.0d0))
        u = (r - transpose(r)) / (2.0d0 * sin(angle)) ! this is skew symmetric
        axis = [u(3,2), u(1,3), u(2,1)]
        axis = axis / norm2(axis)   ! normalize to a unit vector
    end subroutine

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
        !! In vector form this is
        !! $$ \boldsymbol{a}_p = \boldsymbol{a} +
        !! \boldsymbol{\alpha}\times\boldsymbol{x} +
        !! \boldsymbol{\omega}\times(\boldsymbol{\omega}\times\boldsymbol{x}). $$
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
        !! In vector form this is
        !! $$ \boldsymbol{v}_p = \boldsymbol{v} +
        !! \boldsymbol{\omega}\times\boldsymbol{x}. $$
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

end module