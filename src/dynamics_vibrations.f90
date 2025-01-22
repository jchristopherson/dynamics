module dynamics_vibrations
    use iso_fortran_env
    implicit none
    private
    public :: q_factor
    public :: estimate_bandwidth
    public :: logarithmic_decrement

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
end module