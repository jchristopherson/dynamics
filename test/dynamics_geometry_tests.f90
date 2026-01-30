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
    real(real64) :: t, v(3), r0(3), ans(3), s(3)
    type(line) :: ln

    ! Initialization
    call random_number(t)
    call random_number(v)
    call random_number(r0)
    ln%r0 = r0
    ln%v = v
    ans = r0 + t * v
    s = ln%evaluate(t)

    ! Test
    print *, "ANSWER = ", ans
    print *, "COMPUTED = ", s
    ! if (.not.assert(ans, s)) then
    !     rst = .false.
    !     print "(A)", "TEST FAILED: test_line_eval -1"
    ! end if
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
function test_is_parallel_lines() result(rst)
    logical :: rst
    real(real64) :: p1(3), p2(3), np1(3), np2(3), pt(3)
    type(line) :: ln1, ln2, ln3, ln4

    ! Initialization
    rst = .true.
    call random_number(p1)
    p2 = p1 * 3.0   ! parallel to p1
    call random_number(np1)
    call random_number(np2)
    np2(3) = np2(3) + 1.5d0     ! ensuring np1 and np2 are not parallel ever
    call random_number(pt)
    ln1%r0 = pt
    ln1%v = p1
    ln2%r0 = pt
    ln2%v = p2
    ln3%r0 = pt
    ln3%v = np1
    ln4%r0 = pt
    ln4%v = np2

    ! Test
    if (.not.is_parallel(ln1, ln2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_parallel_lines -1"
    end if
    if (is_parallel(ln3, ln4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_parallel_lines -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_is_parallel_planes() result(rst)
    logical :: rst
    real(real64) :: p1(3), p2(3), np1(3), np2(3), pt(3)
    type(plane) :: pln1, pln2, pln3, pln4

    ! Initialization
    rst = .true.
    call random_number(p1)
    p2 = p1 * 3.0   ! parallel to p1
    call random_number(np1)
    call random_number(np2)
    np2(3) = np2(3) + 1.5d0     ! ensuring np1 and np2 are not parallel ever
    call random_number(pt)

    pln1 = plane(pt, p1)
    pln2 = plane(pt, p2)
    pln3 = plane(pt, np1)
    pln4 = plane(pt, np2)

    ! Tests
    if (.not.is_parallel(pln1, pln2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_parallel_planes -1"
    end if
    if (is_parallel(pln3, pln4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_parallel_planes -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_is_point_on_plane() result(rst)
    logical :: rst
    real(real64) :: pt0(3), pt1(3), pt2(3), nrm(3)
    type(plane) :: pln

    ! Initialization
    call random_number(pt0)
    call random_number(pt1)
    call random_number(pt2)
    call random_number(nrm)
    pln = plane(pt0, nrm)
    
    ! Define a point on the plane
    pt1(3) = -(pln%a * pt1(1) + pln%b * pt1(2) + pln%d) / pln%c
    
    ! Define a point off the plane
    pt2 = pt1 + 0.5d0

    ! Test - On Plane
    if (.not.is_point_on_plane(pt1, pln)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_point_on_plane -1"
    end if

    ! Test - Off Plane
    if (is_point_on_plane(pt2, pln)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_point_on_plane -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_point_to_line_distance() result(rst)
    logical :: rst
    real(real64) :: pt1(3), pt2(3), pt(3), t, dist, d
    real(real64), parameter :: t_ans = 0.5d0
    type(line) :: ln

    ! Initialization
    rst = .true.
    call random_number(dist)
    pt1 = [0.0d0, 0.0d0, 0.0d0]
    pt2 = [1.0d0, 1.0d0, 0.0d0]
    pt = t_ans * pt2
    pt(3) = dist
    ln = line(pt1, pt2)

    ! Find the nearest point on the line
    t = nearest_point_on_line(pt, ln)
    if (.not.assert(t, t_ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_point_to_line_distance -1"
    end if

    ! Compute the distance
    d = point_to_line_distance(pt, ln)
    if (.not.assert(dist, d)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_point_to_line_distance -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_is_point_on_line() result(rst)
    logical :: rst
    real(real64) :: pt1(3), pt2(3), ptOn(3), ptOff(3), t
    type(line) :: ln

    ! Initialization
    rst = .true.
    call random_number(t)
    call random_number(pt1)
    call random_number(pt2)
    ln = line(pt1, pt2)
    ptOn = ln%evaluate(t)
    ptOff = ptOn + 0.5d0

    ! Test
    if (.not.is_point_on_line(ptOn, ln)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_point_on_line -1"
    end if

    if (is_point_on_line(ptOff, ln)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_is_point_on_line -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_point_to_plane_distance() result(rst)
    logical :: rst
    real(real64) :: pt(3), pt0(3), nrm(3), dist, d
    type(plane) :: pln

    ! Initialization
    rst = .true.
    call random_number(pt0)
    call random_number(nrm)
    nrm = nrm / norm2(nrm)
    call random_number(dist)
    pt = pt0 + dist * nrm
    pln = plane(pt0, nrm)

    ! Test
    d = point_to_plane_distance(pt, pln)
    if (.not.assert(d, dist)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_point_to_plane_distance -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_vector_plane_projection() result(rst)
    logical :: rst
    real(real64) :: v(3), pt(3), nrm(3), prj(3), nrmans(3)
    type(plane) :: pln

    ! Initialization
    rst = .true.
    pt = [0.0d0, 0.0d0, 0.0d0]
    nrm = [1.0d0, 1.0d0, 1.0d0]
    pln = plane(pt, nrm)
    v = [1.0d0, 0.0d0, 0.0d0]

    prj = vector_plane_projection(v, pln)
end function

! ------------------------------------------------------------------------------
function test_fit_line_to_many_points() result(rst)
    logical :: rst
    real(real64) :: pts(3,3), v(3)
    type(line) :: ln

    ! Initialization
    rst = .true.
    pts = reshape([&
        0.0d0, 0.5d0, 1.0d0, &
        0.0d0, 0.5d0, 1.0d0, &
        0.0d0, 0.5d0, 1.0d0], [3, 3])
    v = pts(3,:) - pts(1,:)
    v = v / norm2(v)
    ln = line(pts)

    ! Test
    if (.not.assert(v, ln%v / norm2(ln%v))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_fit_line_to_many_points -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_fit_plane_to_many_points() result(rst)
    logical :: rst
    real(real64) :: pts(4, 3), nrm(3), pn(3)
    type(plane) :: pln

    ! Initialization
    rst = .true.
    pts = reshape([ &
        0.0d0, 1.0d0, 1.0d0, 0.0d0, &
        0.0d0, 1.0d0, 1.0d0, 0.0d0, &
        0.0d0, 0.0d0, 1.0d0, 1.0d0], [4, 3])
    nrm = [1.0d0, -1.0d0, 0.0d0]
    nrm = nrm / norm2(nrm)
    pln = plane(pts)

    ! Test
    pn = plane_normal(pln)
    if (.not.assert(nrm, pn) .and. .not.assert(nrm, -pn)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_fit_plane_to_many_points -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_flip_plane_normal() result(rst)
    logical :: rst
    real(real64) :: pt(3), nrm(3), n(3)
    type(plane) :: pln

    ! Initialization
    rst = .true.
    call random_number(pt)
    call random_number(nrm)
    nrm = nrm / norm2(nrm)
    pln = plane(pt, nrm)
    call pln%flip_normal()

    ! Test
    n = plane_normal(pln)
    if (.not.assert(-nrm, n)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_flip_plane_normal -1"
    end if
end function

! ------------------------------------------------------------------------------
end module