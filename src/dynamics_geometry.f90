module dynamics_geometry
    use iso_fortran_env
    use dynamics_helper
    implicit none
    private
    public :: plane
    public :: plane_normal
    public :: line
    public :: assignment(=)

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
    end interface

    interface assignment(=)
        module procedure :: plane_assign
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

        rst%a = nrm(1)
        rst%b = nrm(2)
        rst%c = nrm(3)
        rst%d = -dot_product(pt, nrm)
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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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

! ******************************************************************************
! GEOMETRY CALCULATIONS
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module