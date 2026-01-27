module dynamics_maps
    use iso_fortran_env
    use dynamics_geometry
    implicit none
    private
    public :: poincare_map

contains
! ------------------------------------------------------------------------------
    pure function poincare_map(x, y, z, pln) result(rst)
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
        real(real64), allocatable, dimension(:,:) :: rst
            !! An N-by-3 matrix containing the x, y, and z coordinates of each
            !! of the N intersection points in the first, second, and third
            !! columns respectively.

        ! Local Variables
        integer(int32) :: i, j, n
        real(real64) :: t, pt1(3), pt2(3), v(3), pt(3)
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
                j = j + 1
                pt = pt1 + t * v
                buffer(j,:) = pt
            end if

            ! Update the next point
            pt1 = pt2
        end do
        rst = buffer(1:j,:)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module