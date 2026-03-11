module dynamics_linkage_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_linkage
    use dynamics_helper
    use dynamics_rotation
    implicit none
contains
! ------------------------------------------------------------------------------
! Parallel joint axes
function test_define_link_csys_1() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [0.0d0, 0.0d0, 1.0d0]
    rim1 = [0.0d0, 0.0d0, 0.0d0]
    ri = [1.0d0, 1.0d0, 0.0d0]
    u = ri - rim1
    ians = u / norm2(u)
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -4"
    end if
end function

! ------------------------------------------------------------------------------
! Skew, non-intersecting joint axes
function test_define_link_csys_2() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [1.0d0, 0.0d0, 0.0d0]
    rim1 = [0.0d0, 0.0d0, 0.0d0]
    ri = [1.0d0, 1.0d0, 0.0d0]
    u = [0.0d0, 1.0d0, 0.0d0]
    ians = [0.0d0, 1.0d0, 0.0d0]
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -4"
    end if
end function

! ------------------------------------------------------------------------------
! Intersecting joint axes
function test_define_link_csys_3() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [0.0d0, 1.0d0, 0.0d0]
    rim1 = [1.0d0, 1.0d0, 1.0d0]
    ri = [1.0d0, 1.0d0, -1.0d0]
    u = ri
    ians = cross_product(zim1, zi)
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -4"
    end if
end function

! ------------------------------------------------------------------------------
! Colinear joint axes
function test_define_link_csys_4() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [0.0d0, 0.0d0, 1.0d0]
    rim1 = [0.0d0, 0.0d0, 0.0d0]
    ri = [0.0d0, 0.0d0, 1.0d0]
    u = ri
    ians = xim1
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -4"
    end if
end function

! ------------------------------------------------------------------------------
! 2D Serial Manipulator
!
! Jazar's Text
! Example 126
! Page 204
function test_dh_table_1() result(rst)
    logical :: rst
    real(real64) :: l1, l2, l3, t1, t2, t3, R1(3,3), R2(3,3), R3(3,3), R(3,3), &
        i(3), j(3), p1(3), p2(3), p3(3)
    type(coordinate_system) :: csys(4)
    type(dh_table) :: tbl

    ! Parameters
    real(real64), parameter :: ib(3) = [1.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: jb(3) = [0.0d0, 1.0d0, 0.0d0]
    real(real64), parameter :: k(3) = [0.0d0, 0.0d0, 1.0d0]

    ! Initialization
    rst = .true.
    call random_number(l1)
    call random_number(l2)
    call random_number(l3)
    call random_number(t1)
    call random_number(t2)
    call random_number(t3)
    t2 = -t2    ! flip sign on an angle as a sanity check

    ! Define the coordinate systems
    R1 = rotate_z(t1)
    R2 = rotate_z(t2)
    R3 = rotate_z(t3)

    ! CSYS 0
    csys(1) = coordinate_system(ib, jb, k, [0.0d0, 0.0d0, 0.0d0])

    ! CSYS 1
    i = matmul(R1, ib)
    j = matmul(R1, jb)
    p1 = matmul(R1, l1 * ib)
    csys(2) = coordinate_system(i, j, k, p1)

    ! CSYS 2
    R = matmul(R1, R2)
    i = matmul(R, ib)
    j = matmul(R, jb)
    p2 = matmul(R, l2 * ib) + p1
    csys(3) = coordinate_system(i, j, k, p2)

    ! CSYS 3
    R = matmul(R, R3)
    i = matmul(R, ib)
    j = matmul(R, jb)
    p3 = matmul(R, l3 * ib) + p2
    csys(4) = coordinate_system(i, j, k, p3)

    ! Build the table
    tbl = dh_table(csys)

    ! Tests
    if (size(tbl%parameters) /= 3) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1"
        return
    end if

    if (.not.assert(tbl%parameters(1)%link_length, l1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1a"
    end if
    if (.not.assert(tbl%parameters(1)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1b"
    end if
    if (.not.assert(tbl%parameters(1)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1c"
    end if
    if (.not.assert(tbl%parameters(1)%joint_angle, t1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1d"
    end if

    if (.not.assert(tbl%parameters(2)%link_length, l2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2a"
    end if
    if (.not.assert(tbl%parameters(2)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2b"
    end if
    if (.not.assert(tbl%parameters(2)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2c"
    end if
    if (.not.assert(tbl%parameters(2)%joint_angle, t2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2d"
    end if

    if (.not.assert(tbl%parameters(3)%link_length, l3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3a"
    end if
    if (.not.assert(tbl%parameters(3)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3b"
    end if
    if (.not.assert(tbl%parameters(3)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3c"
    end if
    if (.not.assert(tbl%parameters(3)%joint_angle, t3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3d"
    end if
end function

! ------------------------------------------------------------------------------
! 3D PUMA Manipulator
!
! Jazar's Text
! Example 127
! Page 204, 205
function test_dh_table_2() result(rst)
    logical :: rst
    real(real64) :: l1, l2, l3, t1, t2, t3, p1(3), p2(3), p3(3)
    type(coordinate_system) :: csys(4)
    type(dh_table) :: tbl

    ! Parameters
    real(real64), parameter :: ib(3) = [1.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: jb(3) = [0.0d0, 1.0d0, 0.0d0]
    real(real64), parameter :: kb(3) = [0.0d0, 0.0d0, 1.0d0]

    ! Initialization
    rst = .true.
    ! call random_number(l1)
    ! call random_number(l2)
    ! call random_number(l3)
    l1 = 0.5d0
    l2 = 1.5d0
    l3 = 1.0d0
    p1 = [0.0d0, 0.0d0, 0.0d0]
    p2 = [l2, l1, 0.0d0]
    p3 = [0.0d0, 0.0d0, l3] + p2

    ! CSYS 0
    csys(1) = coordinate_system(ib, jb, kb, [0.0d0, 0.0d0, 0.0d0])

    ! CSYS 1
    csys(2) = coordinate_system(ib, kb, -jb, p1)

    ! CSYS 2
    csys(3) = coordinate_system(ib, kb, -jb, p2)

    ! CSYS 3
    csys(4) = coordinate_system(ib, jb, kb, p3)

    ! Build the table
    tbl = dh_table(csys)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module