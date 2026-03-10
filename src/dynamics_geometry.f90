module dynamics_geometry
    use iso_fortran_env
    use dynamics_helper
    use ieee_arithmetic
    use lapack, only : DGESVD, DGELS
    use fstats, only : mean
    implicit none
    private
    public :: plane
    public :: plane_normal
    public :: line
    public :: plucker_line
    public :: assignment(=)
    public :: is_parallel
    public :: is_point_on_plane
    public :: is_point_on_line
    public :: nearest_point_on_line
    public :: point_to_line_distance
    public :: point_to_plane_distance
    public :: vector_plane_projection
    public :: point_plane_projection
    public :: matmul
    public :: line_from_point_and_vector
    public :: line_common_normal
    public :: do_lines_intersect

    type :: plane
        !! Defines a plane as \( a x + b y + c z + d = 0 \).
        real(real64) :: a
            !! The x-component of the plane normal vector.
        real(real64) :: b
            !! The y-component of the plane normal vector.
        real(real64) :: c
            !! The z-component of the plane normal vector.
        real(real64) :: d
            !! The offset from the origin.
    contains
        procedure, public :: flip_normal => plane_flip_normal
    end type

    interface plane
        module procedure :: plane_from_3pts
        module procedure :: plane_from_point_and_normal
        module procedure :: plane_from_many_points
    end interface

    type :: line
        !! Defines the parametric form of a line 
        !! \( \vec{r} = \vec{r_o} + t \vec{v} \).
        real(real64) :: r0(3)
            !! The coordinates of the initial point \(\vec{r_o}\).
        real(real64) :: v(3)
            !! The vector defining the orientation of the line.
    contains
        procedure, public :: evaluate => line_eval
    end type

    interface line
        module procedure :: line_from_2pts
        module procedure :: line_from_2_planes
        module procedure :: line_from_many_points
    end interface

    interface assignment(=)
        module procedure :: plane_assign
        module procedure :: line_assign
        module procedure :: pl_assign
        module procedure :: pl_assign_line
    end interface

    interface is_parallel
        module procedure :: is_parallel_vectors
        module procedure :: is_parallel_lines
        module procedure :: is_parallel_planes
    end interface

    type :: plucker_line
        !! Defines a line in 3D Euclidean space using Plücker coordinates.
        real(real64), public :: v(6)
            !! The 6-element array containing the Plücker coordinates.  The
            !! first 3 elements contain the unit vector and the last 3 elements
            !! contain the moment vector.
    contains
        procedure, public :: u => pl_u
        procedure, public :: m => pl_m
        procedure, public :: to_array => pl_to_array
    end type

    interface plucker_line
        module procedure :: pl_from_2pts
        module procedure :: pl_from_line
        module procedure :: pl_from_2_planes
        module procedure :: pl_from_array
    end interface

    interface matmul
        module procedure :: pl_matmul
    end interface

contains
! ******************************************************************************
! PLANE ROUTINES
! ------------------------------------------------------------------------------
    pure function plane_normal(pln) result(rst)
        !! Returns the normal vector of a plane.
        type(plane), intent(in) :: pln
            !! The plane.
        real(real64) :: rst(3)
            !! The normal vector.

        rst = [pln%a, pln%b, pln%c]
        rst = rst / norm2(rst)
    end function

! ------------------------------------------------------------------------------
    pure function plane_from_3pts(pt1, pt2, pt3) result(rst)
        !! Constructs a plane from 3 points.  The 3 points must not be colinear.
        real(real64), intent(in) :: pt1(3)
            !! The first point.
        real(real64), intent(in) :: pt2(3)
            !! The second point.
        real(real64), intent(in) :: pt3(3)
            !! The third point.
        type(plane) :: rst
            !! The resulting plane.

        real(real64) :: ab(3), ac(3), nrm(3)
        ab = pt2 - pt1
        ac = pt3 - pt1
        nrm = cross_product(ab, ac)
        nrm = nrm / norm2(nrm)
        rst = plane_from_point_and_normal(pt1, nrm)
    end function

! ------------------------------------------------------------------------------
    pure function plane_from_point_and_normal(pt, nrm) result(rst)
        !! Constructs a plane from a point which lies on the plane, and a 
        !! unit vector normal to the plane.
        real(real64), intent(in) :: pt(3)
            !! The point that lies on the plane.
        real(real64), intent(in) :: nrm(3)
            !! The normal unit vector.
        type(plane) :: rst
            !! The resulting plane.

        real(real64) :: nmag
        nmag = norm2(nrm)
        rst%a = nrm(1) / nmag
        rst%b = nrm(2) / nmag
        rst%c = nrm(3) / nmag
        rst%d = -rst%a * pt(1) - rst%b * pt(2) - rst%c * pt(3)
    end function

! ------------------------------------------------------------------------------
    pure function plane_from_many_points(pts) result(rst)
        !! Constructs the plane that best fits a cloud of points in a 
        !! least-squares sense.
        real(real64), intent(in), dimension(:,:) :: pts
            !! The N-by-3 matrix containing the N points to fit.  N must be
            !! at least 3, but is typically much larger.
        type(plane) :: rst
            !! The resulting plane.

        ! Local Variables
        character :: jobu, jobvt
        integer(int32) :: i, m, n, mn, lwork, info
        real(real64), allocatable, dimension(:,:) :: shifted, vt
        real(real64), allocatable, dimension(:) :: s, work
        real(real64) :: temp(1), dummy(1), avg(3), nrm(3), nan

        ! Initialization
        jobu = 'N'
        jobvt = 'S'
        m = size(pts, 1)
        n = size(pts, 2)
        mn = min(m, n)
        nan = ieee_value(nan, IEEE_QUIET_NAN)

        ! Input Checking
        if (m < 3 .or. n /= 3) then
            rst%a = nan
            rst%b = nan
            rst%c = nan
            rst%d = nan
            return
        end if

        ! Workspaec Sizing - only compute V**T
        call DGESVD(jobu, jobvt, m, n, dummy, m, dummy, dummy, m, dummy, mn, &
            temp, -1, info)
        lwork = int(temp(1), int32)
        allocate(work(lwork), vt(mn, n), shifted(m, n), s(mn))

        ! Process
        if (m == 3) then
            ! An exact fit from 3 points
            rst = plane(pts(1,:), pts(2,:), pts(3,:))
        else
            avg(1) = mean(pts(:,1))
            avg(2) = mean(pts(:,2))
            avg(3) = mean(pts(:,3))
            do i = 1, m
                shifted(i,:) = pts(i,:) - avg
            end do
            call DGESVD(jobu, jobvt, m, n, shifted, m, s, dummy, m, vt, mn, &
                work, lwork, info)
            if (info /= 0) then
                rst%a = nan
                rst%b = nan
                rst%c = nan
                rst%d = nan
                return
            end if
            nrm = vt(3,:)
            nrm = nrm / norm2(nrm)
            rst = plane(avg, nrm)
        end if
    end function

! ------------------------------------------------------------------------------
! LINE MEMBER ROUTINES
! ------------------------------------------------------------------------------
    subroutine plane_flip_normal(this)
        !! Flips the normal vector of the plane.
        class(plane), intent(inout) :: this
            !! The plane.

        ! Process
        this%a = -this%a
        this%b = -this%b
        this%c = -this%c
    end subroutine

! ------------------------------------------------------------------------------
! PLANE OPERATORS
! ------------------------------------------------------------------------------
    pure elemental subroutine plane_assign(x, y)
        !! Assigns a plane to another.
        type(plane), intent(out) :: x
            !! The resulting plane.
        type(plane), intent(in) :: y
            !! The source plane

        x%a = y%a
        x%b = y%b
        x%c = y%c
        x%d = y%d
    end subroutine

! ******************************************************************************
! LINE ROUTINES
! ------------------------------------------------------------------------------
    pure function line_from_2pts(pt1, pt2) result(rst)
        !! Constructs a line from two points.
        real(real64), intent(in) :: pt1(3)
            !! The first point.  This point will act as the initial point
            !! along the line such that \( \vec{r_o} = \vec{pt_1}\) in the
            !! equation of the line \( \vec{r} = \vec{r_o} + t \vec{v} \).
        real(real64), intent(in) :: pt2(3)
            !! The second point.
        type(line) :: rst
            !! The resulting line.

        rst%r0 = pt1
        rst%v = pt2 - pt1
    end function

! ------------------------------------------------------------------------------
    pure function line_from_2_planes(p1, p2) result(rst)
        !! Constructs a line from the intersection of two planes.
        class(plane), intent(in) :: p1
            !! The first plane.
        class(plane), intent(in) :: p2
            !! The second plane.
        type(line) :: rst
            !! The resulting line.  NaN's are returned in the event that the
            !! two planes are parallel.

        ! Local Variables
        integer(int32) :: ind
        real(real64) :: n1(3), n2(3), a11, a12, a21, a22, b1, b2, &
            denom, x1, x2

        ! Compute the normal vectors of each plane
        n1 = plane_normal(p1)
        n2 = plane_normal(p2)

        ! Test to see if the planes are parallel
        if (is_parallel(n1, n2)) then
            rst%r0 = ieee_value(0.0d0, IEEE_QUIET_NAN)
            rst%v = ieee_value(0.0d0, IEEE_QUIET_NAN)
            return
        end if

        ! Obtain the direction of the line
        rst%v = cross_product(n1, n2)

        ! Find a point on the line.  The idea is to first locate the index of
        ! the largest magnitue value in the direction vector.  This component
        ! will cross zero at some point, and it is this point we seek to find.
        ind = maxloc(abs(rst%v), 1)
        select case (ind)
        case (1)
            a11 = p1%b
            a21 = p2%b
            a12 = p1%c
            a22 = p2%c
        case (2)
            a11 = p1%a
            a21 = p2%a
            a12 = p1%c
            a22 = p2%c
        case default
            a11 = p1%a
            a21 = p2%a
            a12 = p1%b
            a22 = p2%b
        end select
        b1 = -p1%d
        b2 = -p2%d

        ! Solve the linear system
        denom = a11 * a22 - a12 * a21
        x1 = (a22 * b1 - a12 * b2) / denom
        x2 = (a11 * b2 - a21 * b1) / denom

        ! Find a point along the line to call t = 0
        select case (ind)
        case (1)
            rst%r0(1) = 0.0d0
            rst%r0(2) = x1
            rst%r0(3) = x2
        case (2)
            rst%r0(1) = x1
            rst%r0(2) = 0.0d0
            rst%r0(3) = x2
        case default
            rst%r0(1) = x1
            rst%r0(2) = x2
            rst%r0(3) = 0.0d0
        end select
    end function

! ------------------------------------------------------------------------------
    pure function line_from_many_points(pts) result(rst)
        !! Constructs the line that best fits the supplied set of points in a
        !! least-squares sense.
        real(real64), intent(in), dimension(:,:) :: pts
            !! An N-by-3 matrix where N is at least 2, but typically much
            !! larger.
        type(line) :: rst
            !! The resulting line.

        ! Local Variables
        character :: jobu, jobvt
        integer(int32) :: i, m, n, mn, lwork, info
        real(real64), allocatable, dimension(:,:) :: shifted, vt
        real(real64), allocatable, dimension(:) :: s, work
        real(real64) :: temp(1), dummy(1)

        ! Initialization
        jobu = 'N'
        jobvt = 'S'
        m = size(pts, 1)
        n = size(pts, 2)
        mn = min(m, n)

        ! Input Checking
        if (m < 2 .or. n /= 3) then
            rst%r0 = ieee_value(0.0d0, IEEE_QUIET_NAN)
            rst%v = ieee_value(0.0d0, IEEE_QUIET_NAN)
            return
        end if

        ! Workspace Sizing - only compute V**T
        call DGESVD(jobu, jobvt, m, n, dummy, m, dummy, dummy, m, dummy, mn, &
            temp, -1, info)
        lwork = int(temp(1), int32)
        allocate(work(lwork), vt(mn, n), shifted(m, n), s(mn))

        ! Process
        if (m == 2) then
            rst = line_from_2pts(pts(1,:), pts(2,:))
        else
            rst%r0 = [mean(pts(:,1)), mean(pts(:,2)), mean(pts(:,3))]
            do i = 1, m
                shifted(i,:) = pts(i,:) - rst%r0
            end do
            call DGESVD(jobu, jobvt, m, n, shifted, m, s, dummy, m, vt, mn, &
                work, lwork, info)
            if (info /= 0) then
                rst%r0 = ieee_value(0.0d0, IEEE_QUIET_NAN)
                rst%v = ieee_value(0.0d0, IEEE_QUIET_NAN)
                return
            end if
            rst%v = vt(1,:) / norm2(vt(1,:))
        end if
    end function

! ------------------------------------------------------------------------------
    pure function line_from_point_and_vector(pt, x, nx) result(rst)
        !! Constructs a new line from a point (defines the point where t = 0)
        !! and a direction vector.
        real(real64), intent(in) :: pt(3)
            !! The point at which t = 0 on the line.
        real(real64), intent(in) :: x(3)
            !! The direction vector defining the orientation of the line.
        logical, intent(in), optional :: nx
            !! An optional parameter that defines if x should be normalized to
            !! a unit vector (true), or left as-is (false).  The default is
            !! true such that x is normalized to a unit vector.
        type(line) :: rst
            !! The resulting line.

        ! Local Variables
        logical :: nrm

        ! Initialization
        nrm = .true.
        if (present(nx)) nrm = nx

        ! Process
        rst%r0 = pt
        if (nrm) then
            rst%v = x / norm2(x)
        else
            rst%v = x
        end if
    end function

! ------------------------------------------------------------------------------
! LINE MEMBER ROUTINES
! ------------------------------------------------------------------------------
    pure function line_eval(this, t) result(rst)
        !! Evaluates the equation of the line at the specified parameter.
        class(line), intent(in) :: this
            !! The line.
        real(real64), intent(in) :: t
            !! The parameter.
        real(real64) :: rst(3)
            !! The location along the line defined by the parameter \(t\).

        rst = this%r0 + t * this%v
    end function

! ------------------------------------------------------------------------------
! LINE OPERATORS
! ------------------------------------------------------------------------------
    pure elemental subroutine line_assign(x, y)
        !! Assigns a line to another.
        type(line), intent(out) :: x
            !! The resulting line.
        type(line), intent(in) :: y
            !! The source line

        x%r0 = y%r0
        x%v = y%v
    end subroutine

! ******************************************************************************
! GEOMETRY CALCULATIONS
! ------------------------------------------------------------------------------
    pure function is_parallel_vectors(x, y, tol) result(rst)
        !! Tests to see if two vectors are parallel.
        real(real64), intent(in), dimension(:) :: x
            !! The first vector.
        real(real64), intent(in), dimension(size(x)) :: y
            !! The second vector.
        real(real64), intent(in), optional :: tol
            !! The tolerance to use when testing for parallelism.  The default
            !! tolerance is 10x machine epsilon.
        logical :: rst
            !! Returns true if the vectors are parallel; else, false.

        ! Local Variables
        real(real64) :: t, cp(3)

        ! Initialization
        if (present(tol)) then
            t = tol
        else
            t = 1.0d1 * epsilon(t)
        end if

        ! Process
        cp = cross_product(x, y)
        rst = (abs(cp(1)) <= t .and. abs(cp(2)) <= t .and. abs(cp(3)) <= t)
    end function

! ------------------------------------------------------------------------------
    pure function is_parallel_lines(x, y, tol) result(rst)
        !! Tests to see if two lines are parallel.
        class(line), intent(in) :: x
            !! The first line.
        class(line), intent(in) :: y
            !! The second line.
        real(real64), intent(in), optional :: tol
            !! The tolerance to use when testing for parallelism.  The default
            !! tolerance is 10x machine epsilon.
        logical :: rst
            !! Returns true if the lines are parallel; else, false.

        ! Process
        rst = is_parallel_vectors(x%v, y%v, tol = tol)
    end function

! ------------------------------------------------------------------------------
    pure function is_parallel_planes(x, y, tol) result(rst)
        !! Tests to see if two planes are parallel.
        class(plane), intent(in) :: x
            !! The first plane.
        class(plane), intent(in) :: y
            !! The second plane.
        real(real64), intent(in), optional :: tol
            !! The tolerance to use when testing for parallelism.  The default
            !! tolerance is 10x machine epsilon.
        logical :: rst
            !! Returns true if the planes are parallel; else, false.

        ! Process
        rst = is_parallel_vectors(plane_normal(x), plane_normal(y), tol = tol)
    end function

! ------------------------------------------------------------------------------
    pure function is_point_on_plane(pt, pln, tol) result(rst)
        !! Tests to see if a point lies on a plane.
        real(real64), intent(in) :: pt(3)
            !! The point.
        class(plane), intent(in) :: pln
            !! The plane.
        real(real64), intent(in), optional :: tol
            !! The tolerance to use when testing.  The default
            !! tolerance is 10x machine epsilon.
        logical :: rst
            !! Returns true if the point lies on the plane; else, false.

        ! Local Variables
        real(real64) :: t, s

        ! Initialization
        if (present(tol)) then
            t = tol
        else
            t = 1.0d1 * epsilon(t)
        end if

        ! Process
        s = pln%a * pt(1) + pln%b * pt(2) + pln%c * pt(3) + pln%d
        rst = abs(s) <= t
    end function

! ------------------------------------------------------------------------------
    pure function is_point_on_line(pt, ln, tol) result(rst)
        !! Tests to see if a point lies on a line.
        real(real64), intent(in) :: pt(3)
            !! The point.
        class(line), intent(in) :: ln
            !! The line.
        real(real64), intent(in), optional :: tol
            !! The tolerance to use when testing.  The default
            !! tolerance is 10x machine epsilon.
        logical :: rst
            !! Returns true if the point lies on the line; else, false.

        ! Local Variables
        real(real64) :: t, d

        ! Initialization
        if (present(tol)) then
            t = tol
        else
            t = 1.0d1 * epsilon(t)
        end if

        ! Process
        d = point_to_line_distance(pt, ln)
        rst = d <= t  ! d is always positive
    end function

! ------------------------------------------------------------------------------
    pure function nearest_point_on_line(pt, ln) result(rst)
        !! Gets the line parameter for the point on the line nearest the 
        !! specified point.
        real(real64), intent(in) :: pt(3)
            !! The point.
        class(line), intent(in) :: ln
            !! The line.
        real(real64) :: rst
            !! The line parameteric variable \(t\) defining the location of
            !! the point nearest along the line.

        ! References: 
        ! https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        rst = -dot_product(ln%r0 - pt, ln%v) / (norm2(ln%v)**2)
    end function

! ------------------------------------------------------------------------------
    pure function point_to_line_distance(pt, ln) result(rst)
        !! Computes the shortest distance between a point and a line.
        real(real64), intent(in) :: pt(3)
            !! The point.
        class(line), intent(in) :: ln
            !! The line.
        real(real64) :: rst
            !! The shortest distance between the point and line.

        ! Local Variables
        real(real64) :: t, d(3)

        ! Process
        t = nearest_point_on_line(pt, ln)
        d = pt - ln%evaluate(t)
        rst = norm2(d)
    end function

! ------------------------------------------------------------------------------
    pure function point_to_plane_distance(pt, pln) result(rst)
        !! Computes the shortest distance between a point and a plane.
        real(real64), intent(in) :: pt(3)
            !! The point.
        class(plane), intent(in) :: pln
            !! The plane.
        real(real64) :: rst
            !! The shortest distance between the point and plane.

        ! Process
        rst = abs(pln%a * pt(1) + pln%b * pt(2) + pln%c * pt(3) + pln%d) / &
            norm2([pln%a, pln%b, pln%c])
    end function

! ------------------------------------------------------------------------------
    pure function vector_plane_projection(x, pln) result(rst)
        !! Projects a vector onto a plane.
        real(real64), intent(in) :: x(3)
            !! The vector to project.
        class(plane), intent(in) :: pln
            !! The plane onto which to project the vector.
        real(real64) :: rst(3)
            !! The projected vector.

        ! Local Variables
        real(real64) :: nrm(3), y(3)

        ! Process
        nrm = plane_normal(pln)
        y = vector_projection(x, nrm)
        rst = x - y
    end function

! ------------------------------------------------------------------------------
    pure function point_plane_projection(pt, pln) result(rst)
        !! Projects a point onto a plane.
        real(real64), intent(in) :: pt(3)
            !! The point.
        class(plane), intent(in) :: pln
            !! The plane onto which to project the point.
        real(real64) :: rst(3)
            !! The projected point.

        ! Local Variables
        real(real64) :: t, n(3)

        ! Process
        n = [pln%a, pln%b, pln%c]
        t = -(pln%c * pt(3) + pln%b * pt(2) + pln%a * pt(1) + pln%d) / &
            (norm2(n)**2)
        rst = pt + t * n
    end function

! ******************************************************************************
! PLUCKER_LINE
! ------------------------------------------------------------------------------
    pure function pl_from_2pts(pt1, pt2) result(rst)
        !! Constructs a new plucker_line from two points.
        real(real64), intent(in) :: pt1(3)
            !! The first point.
        real(real64), intent(in) :: pt2(3)
            !! The second point.
        type(plucker_line) :: rst
            !! The resulting line.

        rst%v(1:3) = pt2 - pt1
        rst%v(1:3) = rst%v(1:3) / norm2(rst%v(1:3))
        rst%v(4:6) = cross_product(pt1, rst%v(1:3))
    end function

! ------------------------------------------------------------------------------
    pure function pl_from_line(ln) result(rst)
        !! Constructs a new plucker_line from a line object.
        class(line), intent(in) :: ln
            !! The line.
        type(plucker_line) :: rst
            !! The equivalent plucker_line.

        rst = pl_from_2pts(ln%evaluate(0.0d0), ln%evaluate(1.0d0))
    end function

! ------------------------------------------------------------------------------
    pure function pl_from_2_planes(p1, p2) result(rst)
        !! Constructs a new plucker_line from the intersection of two planes.
        class(plane), intent(in) :: p1
            !! The first plane.
        class(plane), intent(in) :: p2
            !! The second plane.
        type(plucker_line) :: rst
            !! The resulting line.  NaN's are returned in the event that the
            !! two planes are parallel.

        rst = pl_from_line(line(p1, p2))
    end function

! ------------------------------------------------------------------------------
    pure function pl_from_array(x, nrm) result(rst)
        !! Constructs a new plucker_line from the supplied array.
        real(real64), intent(in) :: x(6)
            !! A 6-element array containing the Plücker coordinates.
        logical, intent(in), optional :: nrm
            !! An optional input that specifies if the first three coordinates
            !! (the unit vector) should be normalized (true), or left as-is
            !! (false).  The default is true such that the vector is normalized.
        type(plucker_line) :: rst
            !! The resulting line.

        logical :: n
        n = .true.
        if (present(nrm)) n = nrm
        if (n) then
            rst%v(1:3) = x(1:3) / norm2(x(1:3))
            rst%v(4:6) = x(4:6)
        else
            rst%v = x
        end if
    end function

! ------------------------------------------------------------------------------
    pure function pl_matmul(x, y) result(rst)
        !! Overloads the matmul routine to allow for multiplication of the
        !! Plücker line coordinate vector with a matrix.
        real(real64), intent(in), dimension(:,:) :: x
            !! The N-by-6 matrix.
        type(plucker_line), intent(in) :: y
            !! The plucker_line object.
        real(real64), allocatable, dimension(:) :: rst
            !! The resulting N-element array.

        rst = matmul(x, y%v)
    end function

! ------------------------------------------------------------------------------
! PLUCKER_LINE OPERATORS
! ------------------------------------------------------------------------------
    pure elemental subroutine pl_assign(x, y)
        !! Assigns a plucker_line to another.
        type(plucker_line), intent(out) :: x
            !! The resulting plucker_line.
        type(plucker_line), intent(in) :: y
            !! The source plucker_line.

        x%v = y%v
    end subroutine

! ------------------------------------------------------------------------------
    pure elemental subroutine pl_assign_line(x, y)
        !! Assigns a line to a plucker_line.
        type(plucker_line), intent(out) :: x
            !! The resulting plucker_line.
        type(line), intent(in) :: y
            !! The source line.

        x = plucker_line(y)
    end subroutine

! ------------------------------------------------------------------------------
! PLUCKER_LINE MEMBERS
! ------------------------------------------------------------------------------
    pure function pl_u(this) result(rst)
        !! The unit vector representing the orientation of the line.
        class(plucker_line), intent(in) :: this
            !! The plucker_line object.
        real(real64) :: rst(3)
            !! The unit vector.
        rst = this%v(1:3) 
    end function

! ------------------------------------------------------------------------------
    pure function pl_m(this) result(rst)
        !! The line moment vector.
        class(plucker_line), intent(in) :: this
            !! The plucker_line object.
        real(real64) :: rst(3)
            !! The moment vector.
        rst = this%v(4:6)
    end function

! ------------------------------------------------------------------------------
    pure function pl_to_array(this) result(rst)
        !! Returns the plucker_line as a 6-element array of the form [u, m].
        class(plucker_line), intent(in) :: this
            !! The plucker_line object.
        real(real64) :: rst(6)
            !! The resulting array.

        rst = this%v
    end function

! ******************************************************************************
! ADDITIONAL GEOMETRIC CALCULATIONS (ADDED 3/5/2026, JAC)
! ------------------------------------------------------------------------------
    pure function line_common_normal(ln1, ln2) result(rst)
        !! Returns the common normal line between two lines pointing from ln1 to
        !! ln2.  In the event that the two lines are parallel within the 
        !! specified tolerance, there exist an infinite number of common 
        !! normals; therefore, a line will be chosen that runs from ln1 to ln2
        !! with the point at t = 0 coincident with the point at t = 0 on ln1.
        class(line), intent(in) :: ln1
            !! The first line.
        class(line), intent(in) :: ln2
            !! The second line.
        type(line) :: rst
            !! The common normal line.  The distance along this line between
            !! t = 0 and t = 1 defines the length of the common normal 
            !! connecting ln1 to ln2.  In the event that ln1 and ln2 intersect,
            !! the length of this line is zero.

        ! Local Variables
        logical :: coincident
        integer(int32) :: lwork, info
        real(real64) :: s, t, xi(3), p1(3), p2(3), A(3,2), x(3), temp(1)
        real(real64), allocatable, dimension(:) :: work

        ! Initialization
        t = 1.0d1 * epsilon(t)

        ! Compute the cross products of the direction vectors.
        xi = cross_product(ln2%v, ln1%v)

        ! Locate the point on ln2 that is an intersection between the common 
        ! normal and ln2.  This can be accomplished by noting that:
        !
        ! ln1%r0 + t * xi = ln2%ro + s * ln2%v
        !
        ! We need to solve for s and t
        p1 = ln1%evaluate(0.0d0)

        ! Compute the coefficient matrix
        A(:,1) = xi
        A(:,2) = -ln2%v

        ! Compute the right-hand-side
        x = ln2%r0 - ln1%r0

        ! Set up the least-squares solver.  We use DGELS as we have 3 equations
        ! but only 2 unknowns.
        call DGELS('N', 3, 2, 1, A, 3, x, 3, temp, -1, info)
        lwork = int(temp(1), int32)
        allocate(work(lwork))

        ! Solve
        call DGELS('N', 3, 2, 1, A, 3, x, 3, work, lwork, info)
        if (info > 0) then
            ! The system wasn't of full rank.  It seems the common normal
            ! connecting the two axes is likely of zero length.
            coincident = .true.
        end if
        s = x(2)    ! we can use either t or s.  s is associated with ln2
        p2 = ln2%evaluate(s)
        rst = line(p1, p2)

        ! Ensure that rst is of finite length
        if (norm2(rst%evaluate(1.0d0) - rst%evaluate(0.0d0)) <= t) then
            coincident = .true.
        else if (ieee_is_nan(rst%v(1)) .or. ieee_is_nan(rst%v(2)) .or. &
            ieee_is_nan(rst%v(3))) &
        then
            coincident = .true.
        else
            coincident = .false.
        end if

        if (coincident) then
            ! If the lines are coincident we offset them by the specified
            ! zero tolerance along the computed common normal axis
            rst%r0 = p1
            rst%v = 0.0d0
        end if
    end function

! ------------------------------------------------------------------------------
    pure subroutine do_lines_intersect(ln1, ln2, intersect, t1, t2, tol)
        !! Tests to see if two lines intersect.
        class(line), intent(in) :: ln1
            !! The first line.
        class(line), intent(in) :: ln2
            !! The second line.
        logical, intent(out) :: intersect
            !! True if the two lines intersect within the specified tolerance;
            !! else, false if they do not intersect.
        real(real64), intent(out), optional :: t1
            !! The parametric value associate with ln1 defining the intersection
            !! point.
        real(real64), intent(out), optional :: t2
            !! The parametric value associate with ln2 defining the intersection
            !! point.
        real(real64), intent(in), optional :: tol
            !! The intersection tolerance.  If not supplied, the default value
            !! is 10x machine epsilon.

        ! Local Variables
        integer(int32) :: lwork, info
        real(real64) :: tolerance, A(3, 2), x(3), s, t, temp(1), p1(3), p2(3), &
            dp(3), nan
        real(real64), allocatable, dimension(:) :: work

        ! Initialization
        if (present(tol)) then
            tolerance = tol
        else
            tolerance = 1.0d1 * epsilon(tolerance)
        end if
        nan = ieee_value(nan, IEEE_QUIET_NAN)
        A(:,1) = ln2%v
        A(:,2) = -ln1%v
        x = ln1%r0 - ln2%r0

        ! Set up the least-squares solver
        call DGELS('N', 3, 2, 1, A, 3, x, 3, temp, -1, info)
        lwork = int(temp(1), int32)
        allocate(work(lwork))

        ! Solve A {s;t} = X
        call DGELS('N', 3, 2, 1, A, 3, x, 3, work, lwork, info)
        if (info /= 0) then
            intersect = .false.
            s = nan
            t = nan
        else
            s = x(1)
            t = x(2)
            p1 = ln1%evaluate(t)
            p2 = ln2%evaluate(s)
            dp = abs(p2 - p1)
            if (dp(1) > tolerance .or. dp(2) > tolerance .or. &
                dp(3) > tolerance) &
            then
                ! No intersection
                intersect = .false.
                s = nan
                t = nan
            else
                ! We've got an intersection point
                intersect = .true.
            end if
        end if

        if (present(t1)) t1 = t
        if (present(t2)) t2 = s
    end subroutine

! ------------------------------------------------------------------------------
end module