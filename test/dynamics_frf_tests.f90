module dynamics_frf_tests
    use dynamics
    use iso_fortran_env
    use fortran_test_helper
    implicit none

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

pure subroutine example_2nd_order(t, x, dxdt)
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

pure subroutine example_2nd_order_sweep(freq, x, y, dydx)
    ! Arguments
    real(real64), intent(in) :: freq
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
    f = sin(2.0d0 * pi * freq * x)
    dydx(1) = y(2)
    dydx(2) = f - (2.0d0 * z * wn * y(2) + wn**2 * y(1))
end subroutine
! ------------------------------------------------------------------------------
function test_frf_sweep() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: fmax = 6.0d1
    real(real64), parameter :: fmin = 4.0d1
    real(real64), parameter :: tol = 0.05d0
    integer(int32), parameter :: npts = 10
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    procedure(harmonic_ode), pointer :: fcn
    integer(int32) :: i
    real(real64) :: df, freq(npts)
    real(real64), allocatable, dimension(:) :: mag1, mag2, magans1, magans2, &
        ratio1, ratio2, ref
    type(frf) :: sol
    complex(real64), allocatable, dimension(:) :: tf1, tf2, s

    ! Initialization
    rst = .true.
    fcn => example_2nd_order_sweep
    df = (fmax - fmin) / (npts - 1.0d0)
    freq = (/ (df * i + fmin, i = 0, npts - 1) /)

    ! Compute the FRF's
    sol = frequency_sweep(fcn, freq, [0.0d0, 0.0d0])
    mag1 = abs(sol%responses(:,1))
    mag2 = abs(sol%responses(:,2))

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
subroutine modal_frf_forcing_term(freq, f)
    real(real64), intent(in) :: freq
    complex(real64), intent(out), dimension(:) :: f

    complex(real64), parameter :: zero = (0.0d0, 0.0d0)
    complex(real64), parameter :: one = (1.0d0, 0.0d0)

    f = [one, zero, zero]
end subroutine

! Use the example from: https://github.com/jchristopherson/linalg
function test_proportional_damping_frf() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-3
    integer(int32), parameter :: nfreq = 100
    real(real64), parameter :: fmin = 10.0d0
    real(real64), parameter :: fmax = 1.0d2
    real(real64), parameter :: alpha = 0.1d0
    real(real64), parameter :: beta = 2.0d-5
    real(real64), parameter :: ans1(3) = [232.9225d0, 749.6189d0, 923.5669d0]

    ! Define the model parameters
    real(real64), parameter :: m1 = 0.5d0
    real(real64), parameter :: m2 = 2.5d0
    real(real64), parameter :: m3 = 0.75d0
    real(real64), parameter :: k1 = 5.0d6
    real(real64), parameter :: k2 = 10.0d6
    real(real64), parameter :: k3 = 10.0d6
    real(real64), parameter :: k4 = 5.0d6

    ! Local Variables
    integer(int32) :: i
    real(real64) :: df, m(3,3), k(3,3), freq(nfreq)
    real(real64), allocatable, dimension(:) :: modes
    real(real64), allocatable, dimension(:,:) :: modeshapes
    procedure(modal_excite), pointer :: fcn
    type(frf) :: rsp

    ! Initialization
    rst = .true.
    df = (fmax - fmin) / (nfreq - 1.0d0)
    freq = (/ (df * i + fmin, i = 0, nfreq - 1) /)
    fcn => modal_frf_forcing_term

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])

    ! Compute the frequency response functions, along with the modal information
    rsp = frequency_response(m, k, alpha, beta, freq, fcn, modes, modeshapes)

    ! Verify the modal results
    if (.not.assert(ans1, modes, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_proportional_damping_frf -1"
    end if
end function

! ------------------------------------------------------------------------------
! Use the example from: https://github.com/jchristopherson/linalg
function test_modal_response() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-3
    integer(int32), parameter :: nfreq = 100
    real(real64), parameter :: fmin = 10.0d0
    real(real64), parameter :: fmax = 1.0d2
    real(real64), parameter :: alpha = 0.1d0
    real(real64), parameter :: beta = 2.0d-5
    real(real64), parameter :: ans1(3) = [232.9225d0, 749.6189d0, 923.5669d0]
    real(real64), parameter :: vec1(3) = [0.7179d0, 1.0d0, 0.7466d0]
    real(real64), parameter :: vec2(3) = [-0.4192d0, -0.1638d0, 1.0d0]
    real(real64), parameter :: vec3(3) = [1.0d0, -0.1837d0, 0.1791d0]

    ! Define the model parameters
    real(real64), parameter :: m1 = 0.5d0
    real(real64), parameter :: m2 = 2.5d0
    real(real64), parameter :: m3 = 0.75d0
    real(real64), parameter :: k1 = 5.0d6
    real(real64), parameter :: k2 = 10.0d6
    real(real64), parameter :: k3 = 10.0d6
    real(real64), parameter :: k4 = 5.0d6

    ! Local Variables
    real(real64) :: m(3,3), k(3,3)
    real(real64), allocatable, dimension(:) :: modes
    real(real64), allocatable, dimension(:,:) :: modeshapes

    ! Initialization
    rst = .true.

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])

    ! Compute the modal response
    call modal_response(m, k, modes, modeshapes)

    ! Normalize the mode shapes
    call normalize_mode_shapes(modeshapes)

    ! Verify the modal frequencies
    if (.not.assert(ans1, modes, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_modal_response -1"
    end if

    ! Test the mode shapes
    if (.not.assert(vec1, modeshapes(:,1), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_modal_response -2"
    end if
    if (.not.assert(vec2, modeshapes(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_modal_response -3"
    end if
    if (.not.assert(vec3, modeshapes(:,3), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_modal_response -4"
    end if
end function

! ------------------------------------------------------------------------------
end module