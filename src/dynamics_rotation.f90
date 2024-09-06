module dynamics_rotation
    use iso_fortran_env
    implicit none
    private
    public :: rotate_x
    public :: rotate_y
    public :: rotate_z
    
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
end module