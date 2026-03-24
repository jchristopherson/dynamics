! Compare this example to kinematics_example_1.  This example utilizes the
! serial_linkage type to construct an identical linkage but in a more concise
! manner.

program example
    use iso_fortran_env
    use dynamics
    implicit none
    
    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Properties
    real(real64), parameter :: L1 = 1.5d0
    real(real64), parameter :: L2 = 2.0d1
    real(real64), parameter :: L3 = 1.0d1

    ! Local Variables
    type(binary_link) :: links(3)
    type(serial_linkage) :: linkage
    integer(int32) :: i
    real(real64) :: theta(3), T(4, 4), qo(3), q(3)

    ! Link 1 Definition
    links(1) = binary_link(twist = 0.5d0 * pi, jtype = REVOLUTE_JOINT)
    
    ! Link 2 Definition
    links(2) = binary_link(length = L2, offset = L1, jtype = REVOLUTE_JOINT)

    ! Link 3 Definition
    links(3) = binary_link(length = L3, jtype = REVOLUTE_JOINT)

    ! Build the linkage
    linkage = serial_linkage(links)

    ! Define the joint variables
    call random_number(theta)

    ! --------------------
    ! Compute the forward kinematics & display the matrix
    T = linkage%forward_kinematics(theta)
    do i = 1, 4
        print *, T(i,:)
    end do

    ! --------------------
    ! Solve the inverse problem using the end-effector position and orientation
    ! computed by the forward kinematics process as a target for the inverse
    ! calculations

    ! Define an initial guess for the joint variables
    qo = [0.0d0, 0.0d0, 0.0d0]

    ! Solve the model
    q = linkage%inverse_kinematics(qo, T)

    ! Display the solution and compare with the actual
    print "(A)", "COMPUTED (Inverse Model):"
    print *, q
    print "(A)", "ACTUAL:"
    print *, theta
end program