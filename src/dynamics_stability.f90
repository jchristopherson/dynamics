module dynamics_stability
    use iso_fortran_env
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
function determine_local_stability(a, ev) result(rst)
    !! Determines the nature of stability/unstability near the point at which
    !! the dynamics matrix was computed.
    real(real64), intent(in), dimension(:,:) :: a
        !! An N-by-N matrix containing the 'A' matrix, also known as the
        !! dynamics matrix.
    complex(real64), intent(out), optional, dimension(:) :: ev
        !! An optional N-element array that, if supplied, will be filled with 
        !! the eigenvalues of the matrix A.
    integer(int32) :: rst
        !! Describe the output constants

    ! Local Variables
    logical :: hyperbolic
    integer(int32) :: i, n, npositive, nnegative
    real(real64) :: tol, rv
    complex(real64), allocatable, dimension(:) :: vals
    
    ! Initialization
    n = size(a, 1)
    tol = 1.0d1 * epsilon(tol)  ! zero checking tolerance

    ! Input Checking
    if (size(a, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR

    ! Local Memory Allocation
    allocate(vals(n))

    ! Perform the eigen analysis on A
    call eigen(a, vals)

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
        if (size(ev) /= n) error stop DYN_ARRAY_SIZE_ERROR
        ev = vals
    end if
end function

! ------------------------------------------------------------------------------
end module