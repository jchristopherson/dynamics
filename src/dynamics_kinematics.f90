module dynamics_kinematics
    use iso_fortran_env
    use nonlin_core
    use nonlin_least_squares, only : least_squares_solver
    use ferror
    use dynamics_error_handling
    use dynamics_helper
    implicit none
    private
    public :: identity_4
    public :: dh_rotate_x
    public :: dh_rotate_z
    public :: dh_translate_x
    public :: dh_translate_z
    public :: dh_matrix
    public :: dh_forward_kinematics
    public :: solve_inverse_kinematics
    public :: vecfcn
    public :: least_squares_solver
    public :: iteration_behavior
    public :: jacobian_generating_vector
    public :: dh_jacobian
    public :: REVOLUTE_JOINT
    public :: PRISMATIC_JOINT

    interface dh_forward_kinematics
        module procedure :: dh_forward_kinematics_2
        module procedure :: dh_forward_kinematics_3
        module procedure :: dh_forward_kinematics_4
        module procedure :: dh_forward_kinematics_5
        module procedure :: dh_forward_kinematics_6
        module procedure :: dh_forward_kinematics_7
        module procedure :: dh_forward_kinematics_8
        module procedure :: dh_forward_kinematics_array
    end interface

    interface dh_jacobian
        module procedure :: dh_build_jacobian
    end interface

    integer(int32), parameter :: REVOLUTE_JOINT = 0
        !! Defines a revolute joint.
    integer(int32), parameter :: PRISMATIC_JOINT = 1
        !! Defines a prismatic joint.

! ------------------------------------------------------------------------------
    ! PRIVATE VARIABLES - INVERSE KINEMATICS
    procedure(vecfcn), pointer, private :: kinematics_equations
    real(real64), pointer, private, dimension(:) :: kinematic_constraints

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
    pure function dh_forward_kinematics_array(alpha, a, theta, d) result(rst)
        !! Assembles all of the individual link transformation matrices into a 
        !! single transformation matrix locating the end-effector in the parent
        !! coordinate system for the overall mechanism.  The first entry must
        !! be from the first link nearest ground.
        real(real64), intent(in), dimension(:) :: alpha
            !! The link twist angles, in radians.  This angle is the required
            !! rotation of the z(i-1) axis about the link's x-axis to become
            !! parallel with the link's z-axis.
        real(real64), intent(in), dimension(size(alpha)) :: a
            !! The link lengths as measured along the link's x-axis.
        real(real64), intent(in), dimension(size(alpha)) :: theta
            !! The joint angles, in radians.  This angle is the required rotation
            !! of the z(i-1) axis about the z(i-1) axis to become parallel with
            !! the link's x-axis.
        real(real64), intent(in), dimension(size(alpha)) :: d
            !! The joint offsets distance measured as the distance between the
            !! x(i-1) axis and the link's x-axis along the z(i-1) axis.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: Ti(4,4)

        ! Initialization
        n = size(alpha)
        rst = identity_4()

        ! Process
        do i = 1, n
            Ti = dh_matrix(alpha(i), a(i), theta(i), d(i))
            rst = matmul(rst, Ti)
        end do
    end function

! ------------------------------------------------------------------------------
    function solve_inverse_kinematics(mdl, qo, constraints, df, &
        slvr, ib, err) result(rst)
        !! Solves the inverse kinematics problem for a linkage.  An iterative
        !! solution procedure is utilized.
        procedure(vecfcn), intent(in), pointer :: mdl
            !! A routine used to compute the error in the kinematics 
            !! equations based upon the current solution estimate.
        real(real64), intent(in), dimension(:) :: qo
            !! An M-element array containing an initial estimate of the M joint
            !! variables.
        real(real64), intent(in), target, dimension(:) :: constraints
            !! An N-element array containing the target values (constraints) for
            !! each of the N kinematic equations in the model.  N must be at 
            !! least equal to M (the number of joint variables).
        real(real64), intent(out), optional, target, dimension(:) :: df
            !! An optional N-element array that, if supplied, can be used to 
            !! retrieve the residuals of each of the N kinematic equations.
        class(least_squares_solver), intent(inout), optional, target :: slvr
            !! An optional solver that can be used in place of the default
            !! Levenberg-Marquardt solver.
        type(iteration_behavior), intent(out), optional :: ib
            !! An optional output that can be used to gather information on the
            !! solver.
        class(errors), intent(inout), optional, target :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.
        real(real64), allocatable, dimension(:) :: rst
            !! An M-element array containing the computed joint variables.

        ! Local Variables
        integer(int32) :: nvar, neqn, flag
        real(real64), pointer, dimension(:) :: resid
        real(real64), allocatable, target, dimension(:) :: dresid
        type(vecfcn_helper) :: helper
        class(least_squares_solver), pointer :: solver
        type(least_squares_solver), target :: default_solver
        procedure(vecfcn), pointer :: fcn
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        nvar = size(qo)
        neqn = size(constraints)
        if (present(slvr)) then
            solver => slvr
        else
            solver => default_solver
        end if
        kinematics_equations => mdl
        kinematic_constraints => constraints

        ! Input Check
        if (neqn < nvar) then
            call report_constraint_count_error("solve_inverse_kinematics", &
                nvar, neqn, errmgr)
            return
        end if
        if (present(df)) then
            if (size(df) /= neqn) then
                call report_array_size_error("solve_inverse_kinematics", "df", &
                    neqn, size(df), errmgr)
                return
            end if
        end if

        ! Set up the solver
        fcn => inverse_kinematics_solver
        call helper%set_fcn(fcn, neqn, nvar)

        ! Local Memory Allocations
        allocate(rst(nvar), source = qo, stat = flag)
        if (present(df)) then
            resid => df
        else
            if (flag == 0) allocate(dresid(neqn), stat = flag)
            if (flag == 0) resid => dresid
        end if
        if (flag /= 0) then
            call report_memory_error("solve_inverse_kinematics", flag, errmgr)
            return
        end if

        ! Solve the problem
        call solver%solve(helper, rst, resid, ib = ib, err = errmgr)
    end function

! ----------
    subroutine inverse_kinematics_solver(x, f)
        ! Routine called by the inverse kinematics solver
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f

        ! Compute the kinematics equations
        call kinematics_equations(x, f)

        ! Compare with the constraints
        f = f - kinematic_constraints
    end subroutine

! ******************************************************************************
! V1.0.8 ADDITIONS
! JAN. 29, 2025
! ------------------------------------------------------------------------------
pure function jacobian_generating_vector(d, k, R, jtype) result(rst)
    !! Computes a single Jacobian generating vector given the position vector
    !! of the link origin, \(\vec{d}\), and the joint axis unit vector, 
    !! \(\vec{k}\).  Both vectors must be expressed in the base coordinate
    !! frame.
    !!
    !! $$ \vec{c_{i}} = \left( \begin{matrix}
    !! R \left( \hat{k} \times \vec{d_{i-1}} \right) \\
    !! \vec{k_{i-1}} \end{matrix} \right) $$
    !!
    !! The Jacobian matrix is then constructed from the Jacobian generating
    !! vectors as follows.
    !!
    !! $$ J = \left[ \begin{matrix} \vec{c_1} & \vec{c_2} & ... & \vec{c_n}
    !! \end{matrix} \right] $$
    real(real64), intent(in) :: d(3)
        !! The position vector of the end-effector, \(\vec{d}\), relative to the
        !! link coordinate frame given in the base coordinate frame.
    real(real64), intent(in) :: k(3)
        !! The unit vector defining the joint axis, \(\vec{k}\), given in the
        !! base coordinate frame.
    real(real64), intent(in) :: R(3, 3)
        !! The rotation matrix defining the orientation of the link coordinate
        !! frame relative to the base coordinate frame.
    integer(int32), intent(in) :: jtype
        !! The joint type.  Must be either REVOLUTE_JOINT or PRISMATIC_JOINT.
        !! If incorrectly specified, the code defaults to a REVOLUTE_JOINT type.
    real(real64) :: rst(6)
        !! The resulting 6-element Jacobian generating vector.

    ! Parameter
    real(real64), parameter :: zi(3) = [0.0d0, 0.0d0, 1.0d0]

    ! Local Variables
    real(real64) :: kmag, kcrossd, kunit(3)

    ! Ensure k is a unit vector
    kmag = norm2(k)
    kunit = k / kmag

    ! Process
    if (jtype == PRISMATIC_JOINT) then
        rst(1:3) = kunit
        rst(4:6) = 0.0d0
    else
        ! Compute the cross-product term
        rst(1:3) = matmul(R, cross_product(zi, d))

        ! Fill in the remaining components
        rst(4:6) = kunit
    end if
end function

! ------------------------------------------------------------------------------
function dh_build_jacobian(alpha, a, theta, d, jtypes) result(rst)
    !! Builds the Jacobian matrix for a linkage given the Denavit-Hartenberg
    !! parameters.  The first entry in each array must be from the first link
    !! nearest ground.  The Jacobian matrix relates the joint velocities 
    !! \(\dot{\vec{q}}\) to the end-effector velocity \(\dot{\vec{X}}\) by
    !! \(\dot{\vec{X}} = J \dot{\vec{q}}\).
    real(real64), intent(in), dimension(:) :: alpha
        !! The link twist angles, in radians.  This angle is the required
        !! rotation of the z(i-1) axis about the link's x-axis to become
        !! parallel with the link's z-axis.
    real(real64), intent(in), dimension(size(alpha)) :: a
        !! The link lengths as measured along the link's x-axis.
    real(real64), intent(in), dimension(size(alpha)) :: theta
        !! The joint angles, in radians.  This angle is the required rotation
        !! of the z(i-1) axis about the z(i-1) axis to become parallel with
        !! the link's x-axis.
    real(real64), intent(in), dimension(size(alpha)) :: d
        !! The joint offsets distance measured as the distance between the
        !! x(i-1) axis and the link's x-axis along the z(i-1) axis.
    integer(int32), intent(in), dimension(size(alpha)) :: jtypes
        !! The types of each joint.  Must be either REVOLUTE_JOINT or
        !! PRISMATIC_JOINT.  The code defaults to REVOLUTE_JOINT.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting 6-by-N Jacobian matrix where N is the number of joint
        !! variables (i.e. the length of the input arrays).

    ! Parameters
    real(real64), parameter :: zi_1(4) = [0.0d0, 0.0d0, 1.0d0, 0.0d0]

    ! Local Variables
    integer(int32) :: i, j, n
    real(real64) :: Ti(4, 4), T(4, 4), Te(4, 4), temp(4, 4), di_1(3), ki_1(4)

    ! Initialization
    n = size(alpha)
    T = identity_4()
    allocate(rst(6, n))

    ! Process
    do i = 1, n
        ! Compute the orientation vector
        ki_1 = matmul(T, zi_1)

        ! Compute the link transformation matrix
        Ti = dh_matrix(alpha(i), a(i), theta(i), d(i))

        ! Compute the transformation matrix relating the link to the end 
        ! effector.
        Te = Ti
        do j = i + 1, n
            temp = dh_matrix(alpha(j), a(j), theta(j), d(j))
            Te = matmul(Te, temp)
        end do

        ! Compute the end-effector position vector
        di_1 = Te(1:3,4)

        ! Compute the Jacobian generating vector
        rst(:,i) = jacobian_generating_vector(di_1, ki_1(1:3), T(1:3,1:3), &
            jtypes(i))

        ! Update the transformation matrix
        T = matmul(T, Ti)
    end do
end function

! ------------------------------------------------------------------------------
end module