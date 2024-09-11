module linkage
    use iso_fortran_env
    use dynamics

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Properties
    real(real64), parameter :: L1 = 1.5d0
    real(real64), parameter :: L2 = 2.0d1
    real(real64), parameter :: L3 = 1.0d1

    ! Denavit-Hartenberg Parameters
    real(real64) :: a(3) = [0.0d0, L2, L3]
    real(real64) :: alpha(3) = [0.5d0 * pi, 0.0d0, 0.0d0]
    real(real64) :: d(3) = [0.0d0, L1, 0.0d0]

contains
    subroutine kinematics_equations(jointvars, equations)
        ! The kinematics equations.
        real(real64), intent(in), dimension(:) :: jointvars
            ! The joint variables.
        real(real64), intent(out), dimension(:) :: equations
            ! The resulting kinematic equations.

        ! Local Variables
        real(real64) :: T(4, 4)
        
        ! Compute the forward kinematics problem
        T = dh_forward_kinematics(alpha, a, jointvars, d)

        ! Define the equations.
        ! 1. X position of the end effector
        ! 2. Y position of the end effector
        ! 3. Z position of the end effector
        ! 4. Orientation component of the end-effector
        ! 5. Orientation component of the end-effector
        ! 6. Orientation component of the end-effector
        equations(1:3) = T(1:3,4)
        equations(4) = T(1,1)
        equations(5) = T(2,2)
        equations(6) = T(3,1)
    end subroutine
end module

program example
    use iso_fortran_env
    use dynamics
    use linkage
    implicit none

    ! Local Variables
    integer(int32) :: i
    real(real64) :: theta(3), T(4, 4), qo(3), q(3), constraints(6)
    procedure(vecfcn), pointer :: mdl

    ! Define the Denavit-Hartenberg (DH) parameters
    call random_number(theta)   ! Randomly assign theta.  This is the joint variable

    ! Compute the forward kinematics problem
    T = dh_forward_kinematics(alpha, a, theta, d)

    ! Display the matrix
    do i = 1, 4
        print *, T(i,:)
    end do

    ! -------------------------
    ! Solve the inverse problem.  Use the end-effector position and orientation
    ! computed by the forward kinematics process as a target for the inverse
    ! calculations.

    ! First define an initial guess
    qo = [0.0d0, 0.0d0, 0.0d0]

    ! Define the constraints for each kinematic equation
    constraints(1:3) = T(1:3,4)
    constraints(4) = T(1,1)
    constraints(5) = T(2,2)
    constraints(6) = T(3,1)

    ! Solve the model
    mdl => kinematics_equations
    q = solve_inverse_kinematics(mdl, qo, constraints)

    ! Display the solution and compare with the actual
    print "(A)", "COMPUTED (Inverse Model):"
    print *, q
    print "(A)", "ACTUAL:"
    print *, theta
end program