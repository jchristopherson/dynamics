module dynamics_kinematics
    use iso_fortran_env
    use nonlin
    use ferror
    use dynamics_error_handling
    use dynamics_helper
    use dynamics_geometry
    use dynamics_quaternions
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
    public :: coordinate_system
    public :: dh_parameter_set
    public :: dh_table

    interface dh_forward_kinematics
        module procedure :: dh_forward_kinematics_2
        module procedure :: dh_forward_kinematics_3
        module procedure :: dh_forward_kinematics_4
        module procedure :: dh_forward_kinematics_5
        module procedure :: dh_forward_kinematics_6
        module procedure :: dh_forward_kinematics_7
        module procedure :: dh_forward_kinematics_8
        module procedure :: dh_forward_kinematics_array
        module procedure :: dh_forward_kinematics_dh_params
        module procedure :: dh_forward_kinematics_dh_table
    end interface

    interface dh_jacobian
        module procedure :: dh_build_jacobian
    end interface

    integer(int32), parameter :: REVOLUTE_JOINT = 0
        !! Defines a revolute joint.
    integer(int32), parameter :: PRISMATIC_JOINT = 1
        !! Defines a prismatic joint.

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

! ------------------------------------------------------------------------------
    ! PRIVATE VARIABLES - INVERSE KINEMATICS
    type inverse_kinematics_container
        !! A container for passing kinematics information around the inverse
        !! solver.
        procedure(vecfcn), pointer, nopass :: kinematics_equations
            !! A pointer to the kinematics equations.
        real(real64), allocatable, dimension(:) :: kinematic_constraints
            !! A constraints array.
        class(*), pointer :: user_args => null()
            !! User supplied arguments to the kinematics_equations model.
    end type

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
    pure function dh_forward_kinematics_dh_params(x) result(rst)
        !! Assembles all of the individual link transformation matrices into
        !! a single transformation matrix locating the end-effector in the
        !! parent coordinate system for the overall mechanism.
        class(dh_parameter_set), intent(in), dimension(:) :: x
            !! The list of DH parameters.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: Ti(4, 4)

        ! Initialization
        n = size(x)
        rst = identity_4()

        ! Process
        do i = 1, n
            Ti = dh_matrix(x(i)%link_twist, x(i)%link_length, &
                x(i)%joint_angle, x(i)%link_offset)
            rst = matmul(rst, Ti)
        end do
    end function

! ------------------------------------------------------------------------------
    pure function dh_forward_kinematics_dh_table(x) result(rst)
        !! Assembles all of the individual link transformation matrices into
        !! a single transformation matrix locating the end-effector in the
        !! parent coordinate system for the overall mechanism.
        class(dh_table), intent(in) :: x
            !! The DH table.
        real(real64) :: rst(4, 4)
            !! The resulting 4-by-4 transformation matrix.

        ! Process
        rst = dh_forward_kinematics(x%parameters)
    end function

! ------------------------------------------------------------------------------
    function solve_inverse_kinematics(mdl, qo, constraints, df, &
        slvr, ib, args, err) result(rst)
        !! Solves the inverse kinematics problem for a linkage.  An iterative
        !! solution procedure is utilized.
        procedure(vecfcn), intent(in), pointer :: mdl
            !! A routine used to compute the forward kinematics for the linkage
            !! given the current joint variable estimates.
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
        class(*), intent(inout), optional, target :: args
            !! An optional argument that can be used to communicate with mdl.
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
        type(inverse_kinematics_container) :: obj
        
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
        obj%kinematics_equations => mdl
        obj%kinematic_constraints = constraints
        if (present(args)) obj%user_args => args

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
        call solver%solve(helper, rst, resid, ib = ib, args = obj, err = errmgr)
    end function

! ----------
    subroutine inverse_kinematics_solver(x, f, args)
        ! Routine called by the inverse kinematics solver
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(inout), optional :: args

        ! Process
        select type (args)
        class is (inverse_kinematics_container)
            ! Compute the kinematics equations
            if (associated(args%user_args)) then
                call args%kinematics_equations(x, f, args%user_args)
            else
                call args%kinematics_equations(x, f)
            end if

            ! Compare with the constraints
            f = f - args%kinematic_constraints
        end select
    end subroutine

! ******************************************************************************
! V1.0.8 ADDITIONS
! JAN. 29, 2025
! ------------------------------------------------------------------------------
pure function jacobian_generating_vector(d, k, R, jtype) result(rst)
    !! Computes a single Jacobian generating vector given the position vector
    !! of the link origin, \(\vec{d}\), and the joint axis unit vector, 
    !! \(\vec{k}\).
    !!
    !! For a revolute joint:
    !!
    !! $$ \vec{c_{i}} = \left( \begin{matrix}
    !! R \left( \hat{k} \times \vec{d_{i-1}} \right) \\
    !! \vec{k_{i-1}} \end{matrix} \right) $$
    !!
    !! For a prismatic joint:
    !!
    !! $$ \vec{c_{i}} = \left( \begin{matrix} \vec{k_{i-1}} \\ 0 \end{matrix} 
    !! \right) $$
    !!
    !! The Jacobian matrix is then constructed from the Jacobian generating
    !! vectors as follows.
    !!
    !! $$ J = \left[ \begin{matrix} \vec{c_1} & \vec{c_2} & ... & \vec{c_n}
    !! \end{matrix} \right] $$
    real(real64), intent(in) :: d(3)
        !! The position vector of the end-effector, \(\vec{d}\), relative to the
        !! link coordinate frame given in the base coordinate frame.  An easy
        !! way to compute this vector is to extract the first 3 elements of the
        !! 4th column of the transformation matrix: \(T_{i} T_{i+1} ... T_{n}\).
    real(real64), intent(in) :: k(3)
        !! The unit vector defining the joint axis, \(\vec{k}\), given in the
        !! base coordinate frame.  This vector can be computed most easily by
        !! using the transformation matrix: \(T = T_1 T_2 ... T_{i-1}\) and
        !! then computing \(\vec{k_{i-1}} = T \hat{k}\).
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
            rst = compute_vector_angle(-cs1%k, cs1%i, px2)
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

        ! Parameters
        real(real64), parameter :: half_pi = acos(0.0d0)

        ! Local Variables
        type(quaternion) :: q
        real(real64) :: ry(3), xy, xymag, tol
        logical :: parallel

        ! Initialization
        tol = 1.0d1 * epsilon(tol)

        ! Determine the angle
        rst = vector_angle(x, y)

        ! To determine the sign, we can establish a quaternion and perform the
        ! rotation of one vector onto the other. 
        q = quaternion(rst, axis)   ! angle and axis construction
        ry = aimag(q * y * conjg(q))
        xy = dot_product(x, ry)
        xymag = norm2(x) * norm2(ry)
        parallel = abs(abs(xy) - xymag) < tol
        if (.not.parallel) then
            rst = -rst
        end if

        ! If the angle is pi/2, it is possible for the vectors to point
        ! in opposite directions but still pass the parallelism test
        if (abs(abs(rst) - half_pi) < tol .and. xy < 0.0d0) then
            ! The vectors point in opposite directions
            rst = -rst
        end if
    end function

! ------------------------------------------------------------------------------
end module