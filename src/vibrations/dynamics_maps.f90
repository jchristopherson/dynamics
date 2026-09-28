! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
module dynamics_maps
    use iso_fortran_env
    use dynamics_geometry
    use dynamics_error_handling, only : DYN_INVALID_INPUT_ERROR
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
        !! For consecutive samples \(\boldsymbol{p}_1\) and \(\boldsymbol{p}_2\),
        !! the segment is interpolated as
        !! $$ \boldsymbol{p}(t)=\boldsymbol{p}_1+t(\boldsymbol{p}_2-\boldsymbol{p}_1),
        !! \quad 0\leq t\leq1, $$
        !! and the section point satisfies \(a x(t)+b y(t)+c z(t)+d=0\).
        !! A crossing at an exactly sampled point is returned once, provided
        !! the nearest non-section samples on either side lie on opposite
        !! sides of the plane. Tangencies and runs of samples on the plane are
        !! not crossings and are ignored. A final sample on the plane is
        !! returned once if the preceding segment approaches it.
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
        logical :: from_back, keep
        integer(int32) :: i, j, n, s
        real(real64) :: t, pt(3), normal(3), normal_norm, offset, tol, scale
        real(real64), allocatable, dimension(:) :: signed_distance
        real(real64), allocatable, dimension(:,:) :: buffer
        type(plane) :: p
        
        ! Initialization
        n = size(x)
        allocate(buffer(max(0, n - 1), 3))
        if (n == 0) then
            rst = buffer
            return
        end if
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

        normal = [p%a, p%b, p%c]
        normal_norm = norm2(normal)
        if (normal_norm <= tiny(normal_norm)) error stop DYN_INVALID_INPUT_ERROR
        normal = normal / normal_norm
        offset = p%d / normal_norm
        scale = max(1.0d0, abs(offset), maxval(abs(x)), maxval(abs(y)), maxval(abs(z)))
        tol = 1.0d1 * epsilon(1.0d0) * scale
        allocate(signed_distance(n))
        do i = 1, n
            signed_distance(i) = dot_product(normal, [x(i), y(i), z(i)]) + offset
        end do

        ! Process
        j = 0
        do i = 1, n - 1
            keep = .false.
            if (abs(signed_distance(i)) <= tol) then
                ! An isolated sampled hit belongs to the map only when the
                ! trajectory changes sides across that sample.
                if (i > 1 .and. i < n) then
                    if (abs(signed_distance(i-1)) > tol .and. &
                        abs(signed_distance(i+1)) > tol) then
                        if ((signed_distance(i-1) < 0.0d0 .and. &
                            signed_distance(i+1) > 0.0d0) .or. &
                            (signed_distance(i-1) > 0.0d0 .and. &
                            signed_distance(i+1) < 0.0d0)) then
                            from_back = signed_distance(i-1) < 0.0d0
                            keep = accepts_side(s, from_back)
                            pt = [x(i), y(i), z(i)]
                        end if
                    end if
                end if
            else if (abs(signed_distance(i+1)) > tol) then
                ! Strict opposite signs give a unique segment crossing.
                if ((signed_distance(i) < 0.0d0 .and. &
                    signed_distance(i+1) > 0.0d0) .or. &
                    (signed_distance(i) > 0.0d0 .and. &
                    signed_distance(i+1) < 0.0d0)) then
                    t = signed_distance(i) / &
                        (signed_distance(i) - signed_distance(i+1))
                    pt = [x(i), y(i), z(i)] + t * &
                        ([x(i+1), y(i+1), z(i+1)] - [x(i), y(i), z(i)])
                    from_back = signed_distance(i) < 0.0d0
                    keep = accepts_side(s, from_back)
                end if
            else if (i == n - 1 .and. abs(signed_distance(i)) > tol) then
                ! The final sample has no outgoing segment, but its incoming
                ! direction is known and the endpoint is a valid section hit.
                from_back = signed_distance(i) < 0.0d0
                keep = accepts_side(s, from_back)
                pt = [x(n), y(n), z(n)]
            end if
            if (keep) then
                j = j + 1
                buffer(j,:) = pt
            end if
        end do
        rst = buffer(1:j,:)
    end function

! ------------------------------------------------------------------------------
    pure logical function accepts_side(side, from_back) result(rst)
        integer(int32), intent(in) :: side
        logical, intent(in) :: from_back

        rst = side == POINCARE_TWO_SIDED .or. &
            (side == POINCARE_ONE_SIDED_FROM_BACK .and. from_back) .or. &
            (side == POINCARE_ONE_SIDED_FROM_FRONT .and. .not.from_back)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module