! This example analyzes a planar four-bar linkage.  The mechanism contains a
! single closed kinematic loop, so its forward kinematics require the solution
! of the loop-closure constraints rather than a simple accumulation of link
! transformations.

program example
    use iso_fortran_env
    use dynamics
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    integer(int32), parameter :: npts = 21

    ! Model Properties
    real(real64), parameter :: crank = 1.0d0
    real(real64), parameter :: coupler = 3.5d0
    real(real64), parameter :: rocker = 3.0d0
    real(real64), parameter :: ground = 4.0d0

    ! Local Variables
    integer(int32) :: i
    type(link_container) :: links(4)
    type(joint) :: joints(4)
    type(planar_linkage) :: linkage
    real(real64) :: theta, p(3)
    real(real64), allocatable, dimension(:,:) :: jac

    ! Define the links.  Each link carries a joint frame at each end, with the
    ! body frame located at the first of the two.
    allocate(links(1)%item, source = planar_link(ground))
    allocate(links(2)%item, source = planar_link(crank))
    allocate(links(3)%item, source = planar_link(coupler))
    allocate(links(4)%item, source = planar_link(rocker))

    ! Connect the links.  The crank is the driven member.
    joints(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, actuated = .true.)
    joints(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
    joints(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
    joints(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)

    ! Build the mechanism.  The end-effector is the distal end of the coupler.
    linkage = planar_linkage(links, joints, base = 1, effector = 3, &
        tool = translation(coupler))

    ! Display the mobility of the mechanism
    print "(A,I0)", "Number of independent loops: ", linkage%get_loop_count()
    print "(A,I0)", "Number of joint variables: ", linkage%get_variable_count()
    print "(A,I0)", "Number of constraint equations: ", &
        linkage%get_constraint_count()
    print "(A,I0)", "Degrees of freedom: ", linkage%get_degrees_of_freedom()

    ! Establish a starting configuration.  A closed-loop mechanism admits more
    ! than one assembly mode, so the starting estimate selects the branch.
    call linkage%set_configuration([0.0d0, 0.5d0 * pi, -0.5d0 * pi, 0.0d0])

    ! Sweep the crank and trace the coupler point
    print "(A)", new_line('a') // "  CRANK          X          Y      ANGLE"
    do i = 1, npts
        theta = 2.0d0 * pi * real(i - 1, real64) / real(npts - 1, real64)
        p = linkage%end_effector_pose([theta])
        print "(F7.3, 3F11.4)", theta, p(1), p(2), p(3)
    end do

    ! Examine the sensitivity of the coupler point to the crank at mid-stroke
    theta = 0.5d0 * pi
    jac = linkage%jacobian([theta])
    print "(A)", new_line('a') // "Jacobian at a crank angle of 90 degrees:"
    do i = 1, size(jac, 1)
        print *, jac(i,:)
    end do

contains
    ! Constructs a planar link of the requested length carrying a joint frame at
    ! each end.
    function planar_link(length) result(rst)
        real(real64), intent(in) :: length
        type(multi_joint_link) :: rst
        real(real64) :: frames(4, 4, 2)
        frames(:,:,1) = translation(0.0d0)
        frames(:,:,2) = translation(length)
        rst = multi_joint_link(frames)
    end function

    ! A 4-by-4 translation along the x-axis.
    function translation(a) result(rst)
        real(real64), intent(in) :: a
        real(real64) :: rst(4, 4)
        integer(int32) :: i
        rst = 0.0d0
        do i = 1, 4
            rst(i,i) = 1.0d0
        end do
        rst(1,4) = a
    end function
end program
