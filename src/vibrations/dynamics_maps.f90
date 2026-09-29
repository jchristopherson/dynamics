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
    use diffeq, only : ode_container, ode_integrator, runge_kutta_45
    use dynamics_geometry
    use dynamics_error_handling, only : DYN_INVALID_INPUT_ERROR
    implicit none
    private
    public :: POINCARE_TWO_SIDED
    public :: POINCARE_ONE_SIDED_FROM_FRONT
    public :: POINCARE_ONE_SIDED_FROM_BACK
    public :: poincare_map
    public :: poincare_map_progress

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

    interface poincare_map
        module procedure poincare_map_samples
        module procedure poincare_map_ode
    end interface

    abstract interface
        subroutine poincare_coordinates(t, state, coordinates_out)
            !! Converts an ODE solution sample into coordinates for the
            !! Poincare section. This permits derived coordinates such as
            !! sin(phase) in addition to components of the ODE state.
            import real64
            real(real64), intent(in) :: t
                !! The time at which the ODE state was sampled.
            real(real64), intent(in), dimension(:) :: state
                !! The ODE state at t, in the equation's state ordering.
            real(real64), intent(out), dimension(3) :: coordinates_out
                !! The x, y, and z coordinates to intersect with the plane.
        end subroutine

        subroutine poincare_map_progress(completed_samples, total_samples, &
            time, args)
            !! Reports progress after a complete ODE sample chunk is processed.
            import int32, real64
            integer(int32), intent(in) :: completed_samples
                !! Number of uniformly spaced samples completed so far.
            integer(int32), intent(in) :: total_samples
                !! Total number of requested samples.
            real(real64), intent(in) :: time
                !! Time at the end of the completed chunk.
            class(*), intent(inout), optional :: args
                !! Optional user data shared with the ODE callbacks.
        end subroutine
    end interface

contains
! ------------------------------------------------------------------------------
    pure function poincare_map_samples(x, y, z, pln, side) result(rst)
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
    function poincare_map_ode(sys, tspan, iv, sample_count, pln, side, solver, &
        chunk_size, coordinates, args, progress_callback) result(rst)
        !! Computes a Poincare section from uniformly spaced ODE samples while
        !! retaining only one solution chunk and the resulting section points.
        !! Each chunk starts from the preceding chunk's final solution state.
        !! As with poincare_map_samples, a final sample on the plane can be
        !! retained if the preceding segment approaches it.
        class(ode_container), intent(inout) :: sys
            !! The ODE system to integrate. Its equation function must be set.
        real(real64), intent(in), dimension(2) :: tspan
            !! The increasing start and end times of the complete solve.
        real(real64), intent(in), dimension(:) :: iv
            !! The initial value of each ODE state at tspan(1). At least one
            !! state is required, or three when coordinates is not supplied.
        integer(int32), intent(in) :: sample_count
            !! The number of uniformly spaced samples across tspan, including
            !! both endpoints. Must be at least two.
        class(plane), intent(in), optional :: pln
            !! The section plane. Defaults to the x-y plane (z = 0).
        integer(int32), intent(in), optional :: side
            !! The crossing direction: POINCARE_TWO_SIDED (default),
            !! POINCARE_ONE_SIDED_FROM_FRONT, or POINCARE_ONE_SIDED_FROM_BACK.
        class(ode_integrator), intent(inout), optional, target :: solver
            !! The ODE solver to use. Defaults to runge_kutta_45. Its solution
            !! buffer is cleared for each chunk and on return; other solver
            !! settings, including tolerances, are retained.
        integer(int32), intent(in), optional :: chunk_size
            !! Maximum number of sample intervals per solve. Must be positive;
            !! defaults to 1000. The solve requests at most chunk_size + 1
            !! samples, except that a one-interval solve requests a midpoint.
        procedure(poincare_coordinates), optional :: coordinates
            !! Maps each sampled time and ODE state to section coordinates.
            !! Defaults to the first three state components.
        class(*), intent(inout), optional :: args
            !! Optional user data forwarded to each ODE solver call.
        procedure(poincare_map_progress), intent(in), pointer, optional :: &
            progress_callback
            !! Optional notification after each completed sample chunk.
        real(real64), allocatable, dimension(:,:) :: rst
            !! An N-by-3 array of section intersections in x, y, z order.

        logical :: carry_previous
        integer(int32) :: first, last, count, capacity, chunk, i, n, &
            solve_count, sample_index, offset, found
        real(real64) :: dt, scale, tol, distance, normal_norm
        real(real64), allocatable, dimension(:) :: times, state
        real(real64), allocatable, dimension(:,:) :: solution, points, hits, buffer, copy
        real(real64), dimension(3) :: previous, normal
        type(plane) :: section
        type(runge_kutta_45), target :: default_solver
        class(ode_integrator), pointer :: integrator

        ! Validate the sample grid and coordinate mapping before solving.
        if (sample_count < 2 .or. size(iv) < 1) error stop DYN_INVALID_INPUT_ERROR
        if (.not.present(coordinates) .and. size(iv) < 3) error stop DYN_INVALID_INPUT_ERROR
        if (tspan(2) <= tspan(1)) error stop DYN_INVALID_INPUT_ERROR
        if (.not.sys%get_is_ode_defined()) error stop DYN_INVALID_INPUT_ERROR
        chunk = 1000
        if (present(chunk_size)) chunk = chunk_size
        if (chunk < 1) error stop DYN_INVALID_INPUT_ERROR
        if (present(solver)) then
            integrator => solver
        else
            integrator => default_solver
        end if
        section = plane([0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 1.0d0])
        if (present(pln)) section = pln
        normal = [section%a, section%b, section%c]
        normal_norm = norm2(normal)
        if (normal_norm <= tiny(normal_norm)) error stop DYN_INVALID_INPUT_ERROR
        normal = normal / normal_norm

        ! Allocate storage for one solve and a growable buffer of crossings.
        ! Grid indices are global so chunk boundaries use the same sample times.
        dt = (tspan(2) - tspan(1)) / real(sample_count - 1, real64)
        allocate(times(max(3, min(chunk, sample_count - 1) + 1)), &
            state(size(iv)))
        state = iv
        capacity = 256
        allocate(buffer(capacity, 3))
        count = 0
        first = 0
        carry_previous = .false.
        do while (first < sample_count - 1)
            last = first + min(chunk, sample_count - 1 - first)
            n = last - first + 1
            do i = 1, n
                times(i) = tspan(1) + real(first + i - 1, real64) * dt
            end do
            if (last == sample_count - 1) times(n) = tspan(2)

            ! With only two requested times, DIFFEQ returns every accepted
            ! internal step. Request a midpoint to obtain endpoint samples.
            solve_count = n
            if (n == 2) then
                times(3) = times(2)
                times(2) = 0.5d0 * (times(1) + times(3))
                solve_count = 3
            end if

            ! Release the prior solution before solving the next interval.
            ! The last state becomes the next chunk's initial condition.
            call integrator%clear_buffer()
            call integrator%solve(sys, times(:solve_count), state, args)
            solution = integrator%get_solution()
            if (size(solution,1) /= solve_count) error stop DYN_INVALID_INPUT_ERROR
            state = solution(solve_count,2:)

            ! If the shared endpoint is on the plane, prepend its preceding
            ! sample so the new chunk can distinguish a crossing from a
            ! tangency. Otherwise the shared endpoint alone is sufficient.
            offset = 0
            if (carry_previous) offset = 1
            allocate(points(n + offset, 3))
            if (carry_previous) points(1,:) = previous
            do i = 1, n
                sample_index = i
                if (n == 2) sample_index = 2 * i - 1
                if (present(coordinates)) then
                    call coordinates(solution(sample_index,1), &
                        solution(sample_index,2:), points(i+offset,:))
                else
                    points(i+offset,:) = solution(sample_index,2:4)
                end if
            end do
            hits = poincare_map_samples(points(:,1), points(:,2), points(:,3), &
                section, side)

            ! An interior chunk's final sampled hit is provisional until the
            ! next chunk supplies the sample after it. Match the sample map's
            ! plane normalization and tolerance when testing that endpoint.
            scale = max(1.0d0, abs(section%d / normal_norm), &
                maxval(abs(points)))
            tol = 1.0d1 * epsilon(1.0d0) * scale
            distance = dot_product(normal, points(size(points,1),:)) + &
                section%d / normal_norm
            found = size(hits,1)
            if (last < sample_count - 1 .and. found > 0 .and. &
                abs(distance) <= tol) then
                if (all(hits(found,:) == points(size(points,1),:))) found = found - 1
            end if

            ! Grow storage only with the number of section intersections.
            if (count + found > capacity) then
                capacity = max(2 * capacity, count + found)
                allocate(copy(capacity, 3))
                copy(:count,:) = buffer(:count,:)
                call move_alloc(copy, buffer)
            end if
            buffer(count+1:count+found,:) = hits(:found,:)
            count = count + found

            ! Keep one pre-boundary sample only when the next chunk needs it.
            carry_previous = last < sample_count - 1 .and. abs(distance) <= tol
            if (carry_previous) previous = points(size(points,1)-1,:)
            deallocate(points)
            first = last

            ! Update the user on our progress
            if (present(progress_callback)) then
                call progress_callback(last + 1, sample_count, &
                    times(solve_count), args)
            end if
        end do
        rst = buffer(:count,:)
        call integrator%clear_buffer()
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
end module