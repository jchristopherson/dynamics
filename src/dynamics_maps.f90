module dynamics_maps
    use iso_fortran_env
    use dynamics_geometry
    implicit none
    private
    public :: POINCARE_TWO_SIDED
    public :: POINCARE_ONE_SIDED_FROM_FRONT
    public :: POINCARE_ONE_SIDED_FROM_BACK
    public :: poincare_map

    integer(int32), parameter :: POINCARE_TWO_SIDED = 0
        !! A two-sided Poincare section will be computed.  In this section, the
        !! algorithm does not care whether the trajectory approaches the 
        !! sectioning plane from the front or the back of the plane (defined
        !! by the plane normal).  It simply returns any intersection point.
    integer(int32), parameter :: POINCARE_ONE_SIDED_FROM_FRONT = 1
        !! A one-sided Poincare section will be computed where the algorithm
        !! only retains intersection points where the trajectory approaches
        !! the sectioning plane from the front (the side of the plane normal).
    integer(int32), parameter :: POINCARE_ONE_SIDED_FROM_BACK = 2
        !! A one-sided Poincare section will be computed where the algorithm
        !! only retains intersection points where the trajectory approaches
        !! the sectioning plane from the back (the side opposite the plane 
        !! normal).

contains
! ------------------------------------------------------------------------------
    pure function poincare_map(x, y, z, pln, side) result(rst)
        !! Generates a Poincare map by determining the intersections of the
        !! supplied trajectory with the specified plane.
        real(real64), intent(in), dimension(:) :: x
            !! The x-coordinates of the trajectory.
        real(real64), intent(in), dimension(size(x)) :: y
            !! The y-coordinates of the trajectory.
        real(real64), intent(in), dimension(size(x)) :: z
            !! The z-coordinates of the trajectory.
        class(plane), intent(in), optional :: pln
            !! The plane to intersect.  If not supplied, the x-y plane is 
            !! utilized where z = 0.
        integer(int32), intent(in), optional :: side
            !! An integer flag denoting which approach to use when computing
            !! the section.  The acceptable values are as follows.
            !!
            !! - POINCARE_TWO_SIDED (Default): A two-sided Poincare section 
            !! will be computed.  In this section, the algorithm does not care 
            !! whether the trajectory approaches the sectioning plane from the 
            !! front or the back of the plane (defined by the plane normal).  
            !! It simply returns any intersection point.
            !!
            !! - POINCARE_ONE_SIDED_FROM_FRONT: A one-sided Poincare section 
            !! will be computed where the algorithm only retains intersection 
            !! points where the trajectory approaches the sectioning plane from 
            !! the front (the side of the plane normal).
            !!
            !! - POINCARE_ONE_SIDED_FROM_BACK: A one-sided Poincare section 
            !! will be computed where the algorithm only retains intersection 
            !! points where the trajectory approaches the sectioning plane from 
            !! the back (the side opposite the plane normal).
        real(real64), allocatable, dimension(:,:) :: rst
            !! An N-by-3 matrix containing the x, y, and z coordinates of each
            !! of the N intersection points in the first, second, and third
            !! columns respectively.

        ! Local Variables
        logical :: from_back
        integer(int32) :: i, j, n, s
        real(real64) :: t, pt0(3), pt1(3), pt2(3), v(3), pt(3)
        real(real64), allocatable, dimension(:,:) :: buffer
        type(plane) :: p
        
        ! Initialization
        n = size(x)
        allocate(buffer(n, 3))
        if (present(pln)) then
            p = pln
        else
            ! XY Plane (point & normal)
            p = plane([0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 1.0d0])
        end if
        s = POINCARE_TWO_SIDED
        if (present(side)) then
            if (side == POINCARE_ONE_SIDED_FROM_BACK) then
                s = POINCARE_ONE_SIDED_FROM_BACK
            else if (side == POINCARE_ONE_SIDED_FROM_FRONT) then
                s = POINCARE_ONE_SIDED_FROM_FRONT
            end if
        end if

        ! Process
        j = 0
        pt1 = [x(1), y(1), z(1)]
        do i = 2, n
            ! Construct a line between the two points
            pt2 = [x(i), y(i), z(i)]
            v = pt2 - pt1

            ! Determine the intersection with the plane by determining the
            ! parameter t.  If the parameter t lies between 0 and 1 the point
            ! of intersection is between the two points and we can keep; else, 
            ! it's not and we cannot keep
            t = -(p%a * pt1(1) + p%b * pt1(2) + p%c * pt1(3) + p%d) / &
                dot_product(plane_normal(p), v)
            if (t >= 0.0d0 .and. t <= 1.0d0) then
                pt = pt1 + t * v
                if (s == POINCARE_TWO_SIDED) then
                    j = j + 1
                    buffer(j,:) = pt
                else
                    if (j > 2) then
                        from_back = approach_from_behind(p, pt1, pt0)
                    else
                        from_back = approach_from_behind(p, pt1)
                    end if
                    if (from_back .and. s == POINCARE_ONE_SIDED_FROM_BACK) then
                        ! One sided - from back
                        j = j + 1
                        buffer(j,:) = pt
                    else if (.not.from_back .and. s == POINCARE_ONE_SIDED_FROM_FRONT) then
                        ! One sided - from front
                        j = j + 1
                        buffer(j,:) = pt
                    end if
                end if
            end if

            ! Update the next point
            pt0 = pt1
            pt1 = pt2
        end do
        rst = buffer(1:j,:)
    end function

! ------------------------------------------------------------------------------
    pure function approach_from_behind(pln, x1, x2) result(rst)
        !! Determines if the trajectory approaches the sectioning plane from
        !! behind (true) or not (false) given the most recent point in the
        !! trajectory and the point prior.  The point prior.
        class(plane), intent(in) :: pln
            !! The plane.
        real(real64), intent(in) :: x1(3)
            !! The most recent point in the trajectory prior to the trajectory
            !! intersecting the plane.
        real(real64), intent(in), optional :: x2(3)
            !! The point previos to x1 that is used only if x1 lies on the plane
            !! such that direction is indeterminate.
        logical :: rst
            !! Returns true if the trajectory approaches the plane from behind;
            !! else, false if the trajectory approaches the plane from in front.

        ! Local Variables
        real(real64) :: val, tol

        ! Initialization
        tol = 1.0d1 * epsilon(tol)  ! tolerance for a point on the plane

        ! Process
        val = pln%a * x1(1) + pln%b * x1(2) + pln%c * x1(3) + pln%d
        if (abs(val) < tol) then
            ! The point lies on the plane
            if (present(x2)) then
                ! Use the point previous to determine which side
                val = pln%a * x2(1) + pln%b * x2(2) + pln%c * x2(3) + pln%d
                rst = val < 0.0d0
            else
                rst = .false.   ! default to this case in the event x2 is not given
            end if
        else
            rst = val < 0.0d0   ! the point approaches the plane from the side 
                                ! opposite the normal (from behind) if true
        end if
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module