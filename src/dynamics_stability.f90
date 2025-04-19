module dynamics_stability
    use iso_fortran_env
    use diffeq, only : ode_container
    use ferror
    use dynamics_error_handling
    use linalg, only : eigen
    implicit none
    private
    public :: HYPERBOLIC_FIXED_POINT_SINK
    public :: HYPERBOLIC_FIXED_POINT_SOURCE
    public :: HYPERBOLIC_FIXED_POINT_SADDLE
    public :: NONHYPERBOLIC_FIXED_POINT_UNSTABLE
    public :: NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE
    public :: NONHYPERBOLIC_FIXED_POINT_CENTER
    public :: determine_local_stability

    integer(int32), parameter :: HYPERBOLIC_FIXED_POINT_SINK = 100
        !! Describes a hyperbolic fixed point where all of the eigenvalues of
        !! the dynamics matrix have a nonzero real part and all real parts are
        !! negative-valued.  This point is considered stable.
    integer(int32), parameter :: HYPERBOLIC_FIXED_POINT_SOURCE = 101
        !! Describes a hyperbolic fixed point where all of the eigenvalues of
        !! the dynamics matrix have a nonzero real part and the real
        !! part is positive-valued for each.  This point is considered unstable.
    integer(int32), parameter :: HYPERBOLIC_FIXED_POINT_SADDLE = 102
        !! Describes a hyperbolic fixed point where all of the eigenvalues of
        !! the dynamics matrix have a nonzero real part but one or more of the
        !! eigenvalues has a positive-valued real part.
    integer(int32), parameter :: NONHYPERBOLIC_FIXED_POINT_UNSTABLE = 103
        !! Describes a nonhyperbolic fixed point where one or more of the 
        !! eigenvalues of the dynamics matrix have a positive-valued real part.
    integer(int32), parameter :: NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE = 104
        !! Describes a nonhyperbolic fixed point where some of the eigenvalues 
        !! of the dynamics matrix have negative real parts and the remaining
        !! eigenvalues all have zero-valued real parts.
    integer(int32), parameter :: NONHYPERBOLIC_FIXED_POINT_CENTER = 105
        !! Describes a nonhyperbolic fixed point where all of the eigenvalues
        !! of the dynamics matrix are purely imaginary and nonzero.  This point
        !! is considered stable.

contains
! ------------------------------------------------------------------------------
function determine_local_stability(a, ev, err) result(rst)
    !! Determines the nature of stability/unstability near the point at which
    !! the dynamics matrix was computed.
    real(real64), intent(in), dimension(:,:) :: a
        !! An N-by-N matrix containing the 'A' matrix, also known as the
        !! dynamics matrix.
    complex(real64), intent(out), optional, dimension(:) :: ev
        !! An optional N-element array that, if supplied, will be filled with 
        !! the eigenvalues of the matrix A.
    class(errors), intent(inout), optional, target :: err
        !! An error handler object.
    integer(int32) :: rst
        !! Describe the output constants

    ! Local Variables
    logical :: hyperbolic
    integer(int32) :: i, n, flag, npositive, nnegative
    real(real64) :: tol, rv
    real(real64), allocatable, dimension(:,:) :: acpy
    complex(real64), allocatable, dimension(:) :: vals
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(a, 1)
    tol = 1.0d1 * epsilon(tol)  ! zero checking tolerance

    ! Input Checking
    if (size(a, 2) /= n) then
        call report_nonsquare_matrix_error("determine_local_stability", "a", &
            size(a, 1), size(a, 2), errmgr)
        return
    end if

    ! Local Memory Allocation
    allocate(acpy(n, n), source = a, stat = flag)
    if (flag /= 0) go to 10

    allocate(vals(n), stat = flag)
    if (flag /= 0) go to 10

    ! Perform the eigen analysis on A
    call eigen(acpy, vals, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Cycle over each eigenvalue
    hyperbolic = .true.
    npositive = 0
    nnegative = 0
    do i = 1, n
        rv = real(vals(i), real64)
        if (abs(rv) < tol) then
            ! zero-valued real part - must be nonhyperbolic
            hyperbolic = .false.
        else if (rv > 0.0d0) then
            ! positive-valued real part
            npositive = npositive + 1
        else
            ! negative-valued real part
            nnegative = nnegative + 1
        end if
    end do

    ! Characterize the results
    if (hyperbolic) then
        if (nnegative == n) then
            rst = HYPERBOLIC_FIXED_POINT_SINK
        else if (npositive == n) then
            rst = HYPERBOLIC_FIXED_POINT_SOURCE
        else
            rst = HYPERBOLIC_FIXED_POINT_SADDLE
        end if
    else
        if (nnegative == 0 .and. npositive == 0) then
            rst = NONHYPERBOLIC_FIXED_POINT_CENTER
        else if (nnegative > 0 .and. npositive == 0) then
            rst = NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE
        else
            rst = NONHYPERBOLIC_FIXED_POINT_UNSTABLE
        end if
    end if

    ! Optional Outputs
    if (present(ev)) then
        if (size(ev) /= n) then
            call report_array_size_error("determine_local_stability", "ev", &
                n, size(ev), errmgr)
            return
        end if
        ev = vals
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error("determine_local_stability", flag, errmgr)
    return
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module