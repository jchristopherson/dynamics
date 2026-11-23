program example
    ! This example classifies the stability of several equilibrium points by
    ! linearizing each system about the equilibrium and examining the
    ! eigenvalues of the resulting dynamics matrix.  The first four cases are
    ! a pendulum, governed by theta'' = -(g/L) sin(theta) - c theta', linearized
    ! about its bottom (theta = 0) and top (theta = pi) equilibria.  The last
    ! case is a free mass with viscous damping, x'' = -c x', which has a zero
    ! eigenvalue.
    use iso_fortran_env
    use dynamics
    implicit none

    ! Model Parameters
    real(real64), parameter :: g = 9.81d0
    real(real64), parameter :: pend_length = 1.0d0
    real(real64), parameter :: damping = 5.0d-1
    real(real64), parameter :: negative_damping = -5.0d-1
    real(real64), parameter :: mass_damping = 2.0d0
    real(real64), parameter :: w2 = g / pend_length

! ******************************************************************************
! Stability Analysis
    ! Pendulum with damping at the bottom equilibrium.  Expect a sink.
    call report("Pendulum, damped, bottom equilibrium (theta = 0)", &
        reshape([0.0d0, -w2, 1.0d0, -damping], [2, 2]))

    ! Pendulum with damping at the top equilibrium.  Expect a saddle.
    call report("Pendulum, damped, top equilibrium (theta = pi)", &
        reshape([0.0d0, w2, 1.0d0, -damping], [2, 2]))

    ! Pendulum without damping at the bottom equilibrium.  Expect a center.
    call report("Pendulum, undamped, bottom equilibrium (theta = 0)", &
        reshape([0.0d0, -w2, 1.0d0, 0.0d0], [2, 2]))

    ! Pendulum with negative damping (energy injection) at the bottom
    ! equilibrium.  Expect a source.
    call report("Pendulum, negative damping, bottom equilibrium (theta = 0)", &
        reshape([0.0d0, -w2, 1.0d0, -negative_damping], [2, 2]))

    ! Free mass with viscous damping.  Expect a neutrally stable point.
    call report("Free mass, viscous damping (x' = v, v' = -c v)", &
        reshape([0.0d0, 0.0d0, 1.0d0, -mass_damping], [2, 2]))

contains
    subroutine report(label_, a_)
        !! Classifies the equilibrium described by the dynamics matrix and
        !! prints the classification along with the eigenvalues.
        character(len=*), intent(in) :: label_
            !! A description of the system.
        real(real64), intent(in), dimension(:,:) :: a_
            !! The dynamics matrix.

        ! Local Variables
        integer :: i, flag
        complex(real64), dimension(size(a_, 1)) :: ev

        flag = determine_local_stability(a_, ev = ev)

        write(*,'(A)') label_
        do i = 1, size(ev)
            write(*,'(4X,A,I0,A,2ES13.4)') "lambda_", i, " (real, imag): ", &
                real(ev(i), real64), aimag(ev(i))
        end do
        write(*,'(4X,A,A)') "Classification: ", stability_name(flag)
        write(*,'(A)') ""
    end subroutine

    function stability_name(flag) result(name)
        !! Converts a fixed-point classification constant to a description.
        integer, intent(in) :: flag
            !! The classification constant returned by determine_local_stability.
        character(len=:), allocatable :: name
            !! The description of the classification.

        select case (flag)
        case (HYPERBOLIC_FIXED_POINT_SINK)
            name = "Hyperbolic sink (stable)"
        case (HYPERBOLIC_FIXED_POINT_SOURCE)
            name = "Hyperbolic source (unstable)"
        case (HYPERBOLIC_FIXED_POINT_SADDLE)
            name = "Hyperbolic saddle (unstable)"
        case (NONHYPERBOLIC_FIXED_POINT_UNSTABLE)
            name = "Nonhyperbolic (unstable)"
        case (NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE)
            name = "Nonhyperbolic (neutrally stable)"
        case (NONHYPERBOLIC_FIXED_POINT_CENTER)
            name = "Nonhyperbolic center (stable)"
        case default
            name = "Unknown"
        end select
    end function
end program