module dynamics_geometry
    use iso_fortran_env
    use dynamics_helper
    use ieee_arithmetic
    implicit none
    private
    public :: plane
    public :: plane_normal
    public :: line
    public :: assignment(=)
    public :: is_parallel
    public :: is_point_on_plane
    public :: is_point_on_line
    public :: nearest_point_on_line
    public :: point_to_line_distance

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
    end type

    interface plane
        module procedure :: plane_from_3pts
        module procedure :: plane_from_point_and_normal
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
    end interface

    interface assignment(=)
        module procedure :: plane_assign
        module procedure :: line_assign
    end interface

    interface is_parallel
        module procedure :: is_parallel_vectors
        module procedure :: is_parallel_lines
        module procedure :: is_parallel_planes
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
        ind = maxloc(rst%v, 1)
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
        rst = d <= tol  ! d is always positive
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module