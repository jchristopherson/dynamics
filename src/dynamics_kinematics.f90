module dynamics_kinematics
    use iso_fortran_env
    implicit none
    private
    public :: identity_4
    public :: dh_rotate_x
    public :: dh_rotate_z
    public :: dh_translate_x
    public :: dh_translate_z
    public :: dh_matrix
    public :: dh_forward_kinematics

    interface dh_forward_kinematics
        module procedure :: dh_forward_kinematics_2
        module procedure :: dh_forward_kinematics_3
        module procedure :: dh_forward_kinematics_4
        module procedure :: dh_forward_kinematics_5
        module procedure :: dh_forward_kinematics_6
        module procedure :: dh_forward_kinematics_7
        module procedure :: dh_forward_kinematics_8
    end interface

contains
! ------------------------------------------------------------------------------
    pure function identity_4() result(rst)
        !! Computes a 4-by-4 identity matrix.
        real(real64) :: rst(4, 4)
            !! The resulting identity matrix.

        ! Local Variables & Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Process
        rst(1,1) = one
        rst(2,1) = zero
        rst(3,1) = zero
        rst(4,1) = zero

        rst(1,2) = zero
        rst(2,2) = one
        rst(3,2) = zero
        rst(4,2) = zero

        rst(1,3) = zero
        rst(2,3) = zero
        rst(3,3) = one
        rst(4,3) = zero

        rst(1,4) = zero
        rst(2,4) = zero
        rst(3,4) = zero
        rst(4,4) = one
    end function

! ------------------------------------------------------------------------------
    pure function dh_rotate_x(alpha) result(rst)
        !! Computes the Denavit-Hartenberg matrix for a local x-axis rotation.
        real(real64), intent(in) :: alpha
            !! The rotation angle, in radians.
        real(real64) :: rst(4, 4)
            !! The matrix.

        ! Local Variables & Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64) :: cx, sx

        ! Process
        cx = cos(alpha)
        sx = sin(alpha)

        rst(1,1) = one
        rst(2,1) = zero
        rst(3,1) = zero
        rst(4,1) = zero

        rst(1,2) = zero
        rst(2,2) = cx
        rst(3,2) = sx
        rst(4,2) = zero
        
        rst(1,3) = zero
        rst(2,3) = -sx
        rst(3,3) = cx
        rst(4,3) = zero

        rst(1,4) = zero
        rst(2,4) = zero
        rst(3,4) = zero
        rst(4,4) = one
    end function

! ------------------------------------------------------------------------------
    pure function dh_rotate_z(theta) result(rst)
        !! Computes the Denavit-Hartenberg matrix for a local z-axis rotation.
        real(real64), intent(in) :: theta
            !! The rotation angle, in radians.
        real(real64) :: rst(4, 4)
            !! The matrix.

        ! Local Variables & Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0
        real(real64) :: cx, sx

        ! Process
        cx = cos(theta)
        sx = sin(theta)

        rst(1,1) = cx
        rst(2,1) = sx
        rst(3,1) = zero
        rst(4,1) = zero

        rst(1,2) = -sx
        rst(2,2) = cx
        rst(3,2) = zero
        rst(4,2) = zero
        
        rst(1,3) = zero
        rst(2,3) = zero
        rst(3,3) = one
        rst(4,3) = zero

        rst(1,4) = zero
        rst(2,4) = zero
        rst(3,4) = zero
        rst(4,4) = one
    end function
! ------------------------------------------------------------------------------
    pure function dh_translate_x(a) result(rst)
        !! Computes the Denavit-Hartenberg matrix for a local x-axis 
        !! translation.
        real(real64), intent(in) :: a
            !! The translation.
        real(real64) :: rst(4, 4)
            !! The matrix.

        ! Local Variables & Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Process
        rst(1,1) = one
        rst(2,1) = zero
        rst(3,1) = zero
        rst(4,1) = zero

        rst(1,2) = zero
        rst(2,2) = one
        rst(3,2) = zero
        rst(4,2) = zero

        rst(1,3) = zero
        rst(2,3) = zero
        rst(3,3) = one
        rst(4,3) = zero

        rst(1,4) = a
        rst(2,4) = zero
        rst(3,4) = zero
        rst(4,4) = one
    end function

! ------------------------------------------------------------------------------
    pure function dh_translate_z(d) result(rst)
        !! Computes the Denavit-Hartenberg matrix for a local z-axis 
        !! translation.
        real(real64), intent(in) :: d
            !! The translation.
        real(real64) :: rst(4, 4)
            !! The matrix.

        ! Local Variables & Parameters
        real(real64), parameter :: zero = 0.0d0
        real(real64), parameter :: one = 1.0d0

        ! Process
        rst(1,1) = one
        rst(2,1) = zero
        rst(3,1) = zero
        rst(4,1) = zero

        rst(1,2) = zero
        rst(2,2) = one
        rst(3,2) = zero
        rst(4,2) = zero

        rst(1,3) = zero
        rst(2,3) = zero
        rst(3,3) = one
        rst(4,3) = zero

        rst(1,4) = zero
        rst(2,4) = zero
        rst(3,4) = d
        rst(4,4) = one
    end function

! ------------------------------------------------------------------------------
    pure function dh_matrix(alpha, a, theta, d) result(rst)
        !! Computes the Denavit-Hartenberg transformation matrix for the 
        !! specified DH parameters.
        real(real64), intent(in) :: alpha
            !! The link twist angle, in radians.  This angle is the required
            !! rotation of the z(i-1) axis about the link's x-axis to become
            !! parallel with the link's z-axis.
        real(real64), intent(in) :: a
            !! The link length as measured along the link's x-axis.
        real(real64), intent(in) :: theta
            !! The joint angle, in radians.  This angle is the required rotation
            !! of the z(i-1) axis about the z(i-1) axis to become parallel with
            !! the link's x-axis.
        real(real64), intent(in) :: d
            !! The joint offset distance measured as the distance between the
            !! x(i-1) axis and the link's x-axis along the z(i-1) axis.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Local Variables
        real(real64), dimension(4,4) :: Rx, Dx, Rz, Dz, DxRx, RzDxRx

        ! Compute the matrices
        Rx = dh_rotate_x(alpha)
        Dx = dh_translate_x(a)
        Rz = dh_rotate_z(theta)
        Dz = dh_translate_z(d)

        ! Perform the multiplication
        DxRx = matmul(Dx, Rx)
        RzDxRx = matmul(Rz, DxRx)
        rst = matmul(Dz, RzDxRx)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_2(T1, T2) result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in) :: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Process
        rst = matmul(T1, T2)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_3(T1, T2, T3) result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in) :: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64), intent(in) :: T3(4, 4)
            !! The transformation matrix for the third link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Local Variables
        real(real64) :: T0(4, 4)

        ! Process
        T0 = dh_forward_kinematics_2(T1, T2)
        rst = matmul(T0, T3)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_4(T1, T2, T3, T4) result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in) :: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64), intent(in) :: T3(4, 4)
            !! The transformation matrix for the third link in the linkage.
        real(real64), intent(in) :: T4(4, 4)
            !! The transformation matrix for the fourth link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Local Variables
        real(real64) :: T0(4, 4)

        ! Process
        T0 = dh_forward_kinematics_3(T1, T2, T3)
        rst = matmul(T0, T4)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_5(T1, T2, T3, T4, T5) result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in):: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64), intent(in) :: T3(4, 4)
            !! The transformation matrix for the third link in the linkage.
        real(real64), intent(in) :: T4(4, 4)
            !! The transformation matrix for the fourth link in the linkage.
        real(real64), intent(in) :: T5(4, 4)
            !! The transformation matrix for the fifth link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Local Variables
        real(real64) :: T0(4, 4)

        ! Process
        T0 = dh_forward_kinematics_4(T1, T2, T3, T4)
        rst = matmul(T0, T5)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_6(T1, T2, T3, T4, T5, T6) result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in):: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64), intent(in) :: T3(4, 4)
            !! The transformation matrix for the third link in the linkage.
        real(real64), intent(in) :: T4(4, 4)
            !! The transformation matrix for the fourth link in the linkage.
        real(real64), intent(in) :: T5(4, 4)
            !! The transformation matrix for the fifth link in the linkage.
        real(real64), intent(in) :: T6(4, 4)
            !! The transformation matrix for the sixth link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Local Variables
        real(real64) :: T0(4, 4)

        ! Process
        T0 = dh_forward_kinematics_5(T1, T2, T3, T4, T5)
        rst = matmul(T0, T6)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_7(T1, T2, T3, T4, T5, T6, T7) &
        result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in):: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64), intent(in) :: T3(4, 4)
            !! The transformation matrix for the third link in the linkage.
        real(real64), intent(in) :: T4(4, 4)
            !! The transformation matrix for the fourth link in the linkage.
        real(real64), intent(in) :: T5(4, 4)
            !! The transformation matrix for the fifth link in the linkage.
        real(real64), intent(in) :: T6(4, 4)
            !! The transformation matrix for the sixth link in the linkage.
        real(real64), intent(in) :: T7(4, 4)
            !! The transformation matrix for the seventh link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Local Variables
        real(real64) :: T0(4, 4)

        ! Process
        T0 = dh_forward_kinematics_6(T1, T2, T3, T4, T5, T6)
        rst = matmul(T0, T7)
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_8(T1, T2, T3, T4, T5, T6, T7, T8) &
        result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.
        real(real64), intent(in):: T1(4, 4)
            !! The transformation matrix for the first link nearest ground in
            !! the linkage.
        real(real64), intent(in) :: T2(4, 4)
            !! The transformation matrix for the second link in the linkage.
        real(real64), intent(in) :: T3(4, 4)
            !! The transformation matrix for the third link in the linkage.
        real(real64), intent(in) :: T4(4, 4)
            !! The transformation matrix for the fourth link in the linkage.
        real(real64), intent(in) :: T5(4, 4)
            !! The transformation matrix for the fifth link in the linkage.
        real(real64), intent(in) :: T6(4, 4)
            !! The transformation matrix for the sixth link in the linkage.
        real(real64), intent(in) :: T7(4, 4)
            !! The transformation matrix for the seventh link in the linkage.
        real(real64), intent(in) :: T8(4, 4)
            !! The transformation matrix for the eigth link in the linkage.
        real(real64) :: rst(4, 4)
            !! The resulting transformation matrix.

        ! Local Variables
        real(real64) :: T0(4, 4)

        ! Process
        T0 = dh_forward_kinematics_7(T1, T2, T3, T4, T5, T6, T7)
        rst = matmul(T0, T8)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module