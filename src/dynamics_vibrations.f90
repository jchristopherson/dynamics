module dynamics_vibrations
    use iso_fortran_env
    use peaks
    use ieee_arithmetic
    implicit none
    private
    public :: q_factor
    public :: estimate_bandwidth
    public :: logarithmic_decrement
    public :: damping_from_log_decrement
    public :: find_free_response_properties

contains
! ------------------------------------------------------------------------------
pure function q_factor(zeta) result(rst)
    !! Estimates the Q-factor for a vibratory system.  The Q-factor is computed
    !! \(Q = \frac{1}{2 \zeta}\).
    real(real64), intent(in) :: zeta
        !! The damping ratio.
    real(real64) :: rst
        !! The Q-factor.
    
    ! Process
    rst = 1.0d0 / (2.0d0 * zeta)
end function

! ------------------------------------------------------------------------------
pure function estimate_bandwidth(fn, zeta) result(rst)
    !! Estimates the bandwidth of the resonant mode of a vibratory system.
    !! The bandwidth is the width of the range of frequencies for which the
    !! energy is at least half its peak value and is computed as 
    !! \(\Delta f = \frac{f_n}{Q}\).
    real(real64), intent(in) :: fn
        !! The resonant frequency.  The units are not important; however, 
        !! the units of the output will be the same as the units of this
        !! parameter.
    real(real64), intent(in) :: zeta
        !! The damping ratio.
    real(real64) :: rst
        !! The bandwidth.

    ! Process
    rst = fn / q_factor(zeta)
end function

! ------------------------------------------------------------------------------
pure function logarithmic_decrement(x1, x2, n) result(rst)
    !! Computes the logarithmic decrement given the value of two  successive
    !! peaks in the time history of the free vibratory response of the system.
    !! The logarithmic decrement is calculated as follows.
    !!
    !! $$ \delta = \frac{1}{N} \ln \left( \frac{x(t)}{x(t + N T)} \right) =  
    !! \frac{1}{N} \ln \left( \frac{x_1}{x_2} \right) $$
    real(real64), intent(in) :: x1
        !! The amplitude of the first peak.
    real(real64), intent(in) :: x2
        !! The amplitude of the second peak that occurs N periods after the
        !! first.
    integer(int32), intent(in) :: n
        !! The number of periods of oscillation seperating the two peaks.
    real(real64) :: rst
        !! The logarithmic decrement \(\delta\).

    ! Process
    rst = (1.0d0 / n) * log(x1 / x2)
end function

! ------------------------------------------------------------------------------
pure function damping_from_log_decrement(delta) result(rst)
    !! Computes the damping ratio from the logarithmic decrement \(\delta\).
    !! The damping ratio is related to the logarithmic decrement by the 
    !! following relationship.
    !!
    !! $$ \zeta = \frac{\delta}{\sqrt{4 \pi^2 + \delta^2}} $$
    real(real64), intent(in) :: delta
        !! The logarithmic decrement.
    real(real64) :: rst
        !! The damping ratio.

    ! Process
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    rst = delta / sqrt(4.0d0 * pi**2 + delta**2)
end function

! ------------------------------------------------------------------------------
subroutine find_free_response_properties(t, x, delta, fn, x1, x2, t1, t2, s, n)
    !! Given a free-response time history, this routine attempts to find the 
    !! logarithmic decrement and resonant frequency of a vibratory system. The
    !! logarithmic decrement is estimated by finding successive peaks by
    !! means of peak detection.
    real(real64), intent(in), dimension(:) :: t
        !! An N-element array containing the values in time
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the response sampled at the time points
        !! given in t.
    real(real64), intent(out) :: delta
        !! The logarithmic decrement estimate.  If sufficient peaks cannot be
        !! located, the routine returns NaN.
    real(real64), intent(out) :: fn
        !! The damped resonant frequency in units of Hz, assuming that the
        !! time values are in seconds.  If the time units are not in seconds,
        !! the units will be cycle/unit time with unit time being the units
        !! in which t is supplied.  If sufficient peaks cannot be located, the 
        !! routine returns NaN.
    real(real64), intent(out), optional :: x1
        !! An optional parameter that, if provided, allows for the routine to
        !! return the amplitude of the first peak.  If sufficient peaks cannot 
        !! be located, the routine returns NaN.
    real(real64), intent(out), optional :: x2
        !! An optional parameter that, if provided, allows for the routine to
        !! return the amplitude of the second peak.  If sufficient peaks cannot 
        !! be located, the routine returns NaN.
    real(real64), intent(out), optional :: t1
        !! An optional parameter that, if provided, allows for the routine to
        !! return the time at which the first peak was located.  If sufficient
        !! peaks cannot be located, the routine returns NaN.
    real(real64), intent(out), optional :: t2
        !! An optional parameter that, if provided, allows for the routine to
        !! return the time at which the second peak was located.  If sufficient
        !! peaks cannot be located, the routine returns NaN.
    real(real64), intent(in), optional :: s
        !! An optional input that, if provided, allows for control of the 
        !! sensitivity of the peak detection algorithm.  The default is 0.1%
        !! of the peak-peak amplitude of the signal.
    integer(int32), intent(in), optional :: n
        !! An optional input that, if provided, determines the number of 
        !! periods to allow between peak selection for the logarithmic 
        !! decrement calculation.  The default is 3.

    ! Local Variables
    integer(int32) :: np, i1, i2, j2
    real(real64) :: xmax, xmin, dx, x1p, x2p, t1p, t2p, nan
    integer(int32), allocatable, dimension(:) :: maxind, minind
    real(real64), allocatable, dimension(:) :: maxvals, minvals

    ! Determine a suitable sensitivity to peak detection
    if (present(s)) then
        dx = s
    else
        xmax = maxval(x)
        xmin = minval(x)
        dx = 1.0d-3 * (xmax - xmin)
    end if

    ! Peak Count
    if (present(n)) then
        np = n
    else
        np = 3
    end if

    ! Additional initialization
    nan = ieee_value(nan, IEEE_QUIET_NAN)
    delta = nan
    fn = nan
    t1p = nan
    t2p = nan
    x1p = nan
    x2p = nan

    ! Locate peaks
    call peak_detect(x, dx, maxind, maxvals, minind, minvals)
    if (size(maxind) < 2) then
        ! Return NaN's as we couldn't find enough peaks
        go to 10
    end if
    i1 = maxind(1)
    np = min(np, size(maxind) - 1)
    j2 = np + 1
    i2 = maxind(j2)
    t1p = t(i1)
    t2p = t(i2)
    x1p = x(i1)
    x2p = x(i2)
    delta = logarithmic_decrement(x1p, x2p, np)
    fn = np / (t2p - t1p)

    ! End
10  continue
    if (present(x1)) x1 = x1p
    if (present(x2)) x2 = x2p
    if (present(t1)) t1 = t1p
    if (present(t2)) t2 = t2p
end subroutine

! ------------------------------------------------------------------------------
end module