module dynamics_linkage_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_linkage
    use dynamics_helper
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

! ------------------------------------------------------------------------------
end module