module dynamics_linkage_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_linkage
    use dynamics_helper
    use dynamics_rotation
    use dynamics_kinematics
    implicit none
contains
! ------------------------------------------------------------------------------
    function build_serial_linkage_1() result(rst)
        !! Builds a spherical manipulator (Jazar's text example 214, p 353).
        type(serial_linkage) :: rst

        ! Model Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: l0 = 1.5d0
        real(real64), parameter :: dstatic = 1.25d0

        ! Local Variables
        type(binary_link) :: links(3)
        real(real64) :: J(6,3)

        ! Initialization
        links(1) = binary_link(twist = -0.5d0 * pi, offset = l0)
        links(2) = binary_link(twist = 0.5d0 * pi)
        links(3) = binary_link(offset = dstatic, jtype = PRISMATIC_JOINT)
        
        ! Build the linkage
        rst = serial_linkage(links)
    end function

! ------------------------------------------------------------------------------
    function test_serial_linkage_forward_kinematics() result(rst)
        ! Arguments
        logical :: rst
        
        ! Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        integer(int32) :: i
        class(binary_link), pointer :: lnk
        type(serial_linkage) :: linkage
        real(real64), dimension(3) :: q, d, alpha, a, theta
        real(real64), dimension(4, 4) :: T1, T2a, T2, T3a, Tans, T

        ! Initialization
        rst = .true.
        linkage = build_serial_linkage_1()  ! assuming a 3-element linkage
        call random_number(q)
        do i = 1, 3
            lnk => linkage%get_link(i)
            alpha(i) = lnk%link_twist
            a(i) = lnk%link_length
            if (lnk%joint_type == PRISMATIC_JOINT) then
                theta(i) = lnk%joint_angle
                d(i) = lnk%link_offset + q(i)
            else
                theta(i) = lnk%joint_angle + q(i)
                d(i) = lnk%link_offset
            end if
        end do

        ! Compute the solution
        T1 = dh_matrix(alpha(1), a(1), theta(1), d(1))
        T2a = dh_matrix(alpha(2), a(2), theta(2), d(2))
        T2 = matmul(T1, T2a)
        T3a = dh_matrix(alpha(3), a(3), theta(3), d(3))
        Tans = matmul(T2, T3a)

        ! Evaluate the linkage forward kinematics
        T = linkage%forward_kinematics(q)

        ! Test
        if (.not.assert(T, Tans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_serial_linkage_forward_kinematics -1"
        end if
    end function

! ------------------------------------------------------------------------------
end module