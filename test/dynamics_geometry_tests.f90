module dynamics_geometry_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_geometry
    use dynamics_helper
    implicit none

contains
! ------------------------------------------------------------------------------
function test_line_from_2_points() result(rst)
    logical :: rst
    real(real64) :: pt1(3), pt2(3), v(3)
    type(line) :: ln

    ! Initialization
    rst = .true.
    call random_number(pt1)
    call random_number(pt2)
    v = pt2 - pt1
    ln = line(pt1, pt2)

    ! Tests
    if (.not.assert(pt1, ln%r0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_line_from_2_points -1"
    end if

    if (.not.assert(v, ln%v)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_line_from_2_points -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_line_from_2_planes() result(rst)
    logical :: rst
    real(real64) :: n1(3), n2(3), v(3), pt(3)
    type(plane) :: p1, p2
    type(line) :: ln

    ! Initialization
    rst = .true.
    n1 = [1.0d0, 1.0d0, 0.0d0]
    n2 = [0.0d0, 0.0d0, 1.0d0]
    v = [1.0d0, -1.0d0, 0.0d0]
    pt = [0.0d0, 0.0d0, 0.0d0]

    n1 = n1 / norm2(n1)
    n2 = n2 / norm2(n2)
    v = v / norm2(v)

    p1 = plane(pt, n1)
    p2 = plane(pt, n2)
    ln = line(p1, p2)

    ! Tests
    if (.not.assert(pt, ln%r0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_line_from_2_planes -1"
    end if

    if (.not.assert(v, ln%v)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_line_from_2_planes -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_line_eval() result(rst)
    logical :: rst
    real(real64) :: t, v(3), r0(3)
    type(line) :: ln

    ! Initialization
    call random_number(t)
    call random_number(v)
    call random_number(r0)
    ln%r0 = r0
    ln%v = v

    ! Test
    if (.not.assert(r0 + t * v, ln%evaluate(t))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_line_eval -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_plane_from_3_points() result(rst)
    logical :: rst
    real(real64) :: pt1(3), pt2(3), pt3(3), nrm(3)
    type(plane) :: pln

    ! Initialization
    rst = .true.
    call random_number(pt1)
    call random_number(pt2)
    call random_number(pt3)
    nrm = cross_product(pt2 - pt1, pt3 - pt1)
    nrm = nrm / norm2(nrm)
    pln = plane(pt1, pt2, pt3)

    ! Test
    if (.not.assert(pln%a, nrm(1))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_3_points -1"
    end if
    if (.not.assert(pln%b, nrm(2))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_3_points -2"
    end if
    if (.not.assert(pln%c, nrm(3))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_3_points -3"
    end if
    if (.not.assert(pln%a * pt1(1) + pln%b * pt1(2) + pln%c * pt1(3), -pln%d)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_3_points -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_plane_from_point_and_normal() result(rst)
    logical :: rst
    real(real64) :: pt1(3), nrm(3)
    type(plane) :: pln

    ! Initialization
    rst = .true.
    call random_number(pt1)
    call random_number(nrm)
    nrm = nrm / norm2(nrm)
    pln = plane(pt1, nrm)

    ! Test
    if (.not.assert(pln%a, nrm(1))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_point_and_normal -1"
    end if
    if (.not.assert(pln%b, nrm(2))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_point_and_normal -2"
    end if
    if (.not.assert(pln%c, nrm(3))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_point_and_normal -3"
    end if
    if (.not.assert(pln%a * pt1(1) + pln%b * pt1(2) + pln%c * pt1(3), -pln%d)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_from_point_and_normal -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_plane_normal() result(rst)
    logical :: rst
    real(real64) :: pt1(3), nrm(3)
    type(plane) :: pln

    ! Initialization
    rst = .true.
    call random_number(pt1)
    call random_number(nrm)
    pln = plane(pt1, nrm)

    ! Test
    nrm = nrm / norm2(nrm)
    if (.not.assert(nrm, plane_normal(pln))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_plane_normal -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_is_parallel_vectors() result(rst)
    logical :: rst
    real(real64) :: p1(3), p2(3), np1(3), np2(3)

    ! Initialization
    rst = .true.
    call random_number(p1)
    p2 = p1 * 3.0   ! parallel to p1
    call random_number(np1)
    call random_number(np2)
    np2(3) = np2(3) + 1.5d0     ! ensuring np1 and np2 are not parallel ever

    print *, cross_product(p1, p2)
    if (.not.is_parallel(p1, p2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_parallel_vectors -1"
    end if

    if (is_parallel(np1, np2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_parallel_vectors -2"
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module