module dynamics_frf_tests
    use dynamics
    use iso_fortran_env
    use fortran_test_helper
    implicit none

    type, extends(harmonic_ode_container) :: ode_sweep_object
    contains
        procedure, public :: ode => example_2nd_order_sweep
    end type

contains
! ------------------------------------------------------------------------------
! 2nd Order Test Problem
! x" + 2 * z * wn * x' + wn**2 * x = f(t)
pure function example_2nd_order_forcing(t) result(rst)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64) :: rst

    ! Process
    rst = chirp(t, 1.0d2, 5.0d0, 1.0d0, 1.0d2)
end function

subroutine example_2nd_order(t, x, dxdt)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: dxdt

    ! Model Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1

    ! Local Variables
    real(real64) :: f

    ! Process
    f = example_2nd_order_forcing(t)
    dxdt(1) = x(2)
    dxdt(2) = f - (2.0d0 * z * wn * x(2) + wn**2 * x(1))
end subroutine

subroutine example_2nd_order_sweep(this, x, y, dydx)
    ! Arguments
    class(ode_sweep_object), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx

    ! Model Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1

    ! Local Variables
    real(real64) :: f

    ! Process
    f = sin(2.0d0 * pi * this%excitation_frequency * x)
    dydx(1) = y(2)
    dydx(2) = f - (2.0d0 * z * wn * y(2) + wn**2 * y(1))
end subroutine

! ------------------------------------------------------------------------------
function test_fft_based_frf() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: fmax = 1.0d2
    real(real64), parameter :: tol = 0.05d0
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    type(forced_ode_container) :: mdl
    integer(int32) :: i, ind, start
    real(real64), allocatable, dimension(:) :: freq, mag1, mag2, magans1, &
        magans2, ratio, reference, reference2
    complex(real64), allocatable, dimension(:,:) :: sol
    complex(real64), allocatable, dimension(:) :: tf1, tf2, s

    ! Initialization
    rst = .true.
    mdl%fcn => example_2nd_order
    mdl%forcing_function => example_2nd_order_forcing
    
    ! Compute the FRF's
    call frequency_response(mdl, 5.0d0, [0.0d0, 0.0d0], sol, freq = freq)
    mag1 = abs(sol(:,1))
    mag2 = abs(sol(:,2))

    ! Compute the solution
    s = j * (2.0d0 * pi * freq)
    tf1 = 1.0d0 / (s**2 + 2.0d0 * z * wn * s + wn**2)
    tf2 = s * tf1
    magans1 = abs(tf1)
    magans2 = abs(tf2)

    ! The solution will carry on beyond the desired max frequency; however, is
    ! effectively invalid beyond such frequency.  As a result we must limit the
    ! test range
    ind = 0
    do i = 1, size(mag1)
        if (freq(i) > fmax) then
            ind = i - 1
            exit
        end if
    end do

    ! Test
    ratio = mag1(:ind) / magans1(:ind)
    allocate(reference(size(ratio)), source = 1.0d0)
    if (.not.assert(ratio, reference, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_fft_based_frf -1"
    end if

    start = floor(0.1 * ind)
    ratio = mag2(start:ind) / magans2(start:ind)  ! The non-one start is because for small s-values small errors can lead to large ratio differences
    allocate(reference2(size(ratio)), source = 1.0d0)
    if (.not.assert(ratio, reference2, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_fft_based_frf -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_frf_sweep() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: fmax = 1.0d2
    real(real64), parameter :: fmin = 1.0d1
    real(real64), parameter :: tol = 0.05d0
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    type(ode_sweep_object) :: mdl
    integer(int32) :: i
    real(real64) :: df, freq(npts)
    real(real64), allocatable, dimension(:) :: mag1, mag2, magans1, magans2, &
        ratio1, ratio2, ref
    complex(real64), allocatable, dimension(:,:) :: sol
    complex(real64), allocatable, dimension(:) :: tf1, tf2, s

    ! Initialization
    rst = .true.
    df = (fmax - fmin) / (npts - 1.0d0)
    freq = (/ (df * i + fmin, i = 0, npts - 1) /)

    ! Compute the FRF's
    sol = mdl%frequency_sweep(freq, [0.0d0, 0.0d0])
    mag1 = abs(sol(:,1))
    mag2 = abs(sol(:,2))

    ! Compute the solution
    s = j * (2.0d0 * pi * freq)
    tf1 = 1.0d0 / (s**2 + 2.0d0 * z * wn * s + wn**2)
    tf2 = s * tf1
    magans1 = abs(tf1)
    magans2 = abs(tf2)

    ! Tests
    ratio1 = mag1 / magans1
    allocate(ref(size(ratio1)), source = 1.0d0)
    if (.not.assert(ratio1, ref, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_frf_sweep -1"
    end if

    ratio2 = mag2 / magans2
    if (.not.assert(ratio2, ref, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_frf_sweep -2"
    end if
end function

! ------------------------------------------------------------------------------
end module