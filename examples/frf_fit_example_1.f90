module excitation_fit_example
    use iso_fortran_env
    implicit none

contains
    subroutine modal_frf_forcing_term(freq, f)
        real(real64), intent(in) :: freq
        complex(real64), intent(out), dimension(:) :: f

        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        f = [1.0d3 * one, zero, zero]
    end subroutine
end module

program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use excitation_fit_example
    implicit none

    ! Parameters
    integer(int32), parameter :: nfreq = 1000
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: fmin = 2.0d0 * pi * 10.0d0
    real(real64), parameter :: fmax = 2.0d0 * pi * 1.0d3
    real(real64), parameter :: alpha = 1.0d-3
    real(real64), parameter :: beta = 2.0d-6
    character(len = *), parameter :: tab = achar(9)

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
    real(real64) :: m(3,3), k(3,3)
    type(frf) :: rsp
    procedure(modal_excite), pointer :: fcn
    real(real64), allocatable, dimension(:) :: mdl, freq, amp, mdlamp, phase, &
        mdlphase, natfreqs
    complex(real64), allocatable, dimension(:) :: mdlrsp

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: x1, x2, y1, y2
    class(legend), pointer :: lgnd

    ! Initialization
    fcn => modal_frf_forcing_term

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])

    ! Compute the frequency response functions
    rsp = frequency_response(m, k, alpha, beta, nfreq, fmin, fmax, fcn)

    ! Attempt to fit the modeled frequency response function
    mdl = fit_frf(FRF_FORCE_EXCITATION, 3, rsp%frequency, rsp%responses(:,1))

    ! Evaluate the model
    mdlrsp = evaluate_frf_force_model(mdl, rsp%frequency)

    ! Extract amplitude and phase information
    amp = 2.0d1 * log10(abs(rsp%responses(:,1) / rsp%responses(1,1)))
    phase = 1.8d2 * atan2(aimag(rsp%responses(:,1)), real(rsp%responses(:,1))) / pi

    mdlamp = 2.0d1 * log10(abs(mdlrsp / mdlrsp(1)))
    mdlphase = 1.8d2 * atan2(aimag(mdlrsp), real(mdlrsp)) / pi

    ! Get the frequency values
    freq = rsp%frequency / (2.0d0 * pi)

    ! Print out the model coefficients
    print "(A)", "MODEL COEFFICIENTS:"
    do i = 1, size(mdl) / 3
        print "(AI0AF0.4)", tab // "A", i, ": ", mdl(3 * i - 2)
        print "(AI0AF0.4A)", tab // "wn", i, ": ", &
            mdl(3 * i - 1) / (2.0d0 * pi), " Hz"
        print "(AI0AF0.4)", tab // "z", i, ": ", mdl(3 * i)
    end do

    ! For comparison purposes, compute and display the actual modal 
    ! frequencies of the system.
    call modal_response(m, k, natfreqs)
    natfreqs = natfreqs / (2.0d0 * pi)  ! convert units to Hz
    print "(A)", new_line('a') // "ACTUAL MODAL FREQUENCIES:"
    do i = 1, size(natfreqs)
        print "(AI0AF0.4A)", tab // "Mode ", i, ": ", natfreqs(i), " Hz"
    end do

! ------------------------------------------------------------------------------
! PLOTTING CODE ONLY
! ------------------------------------------------------------------------------
    ! Set up the plots
    call plt%initialize(2, 1)
    call plt1%initialize()
    call plt2%initialize()
    x1 => plt1%get_x_axis()
    y1 => plt1%get_y_axis()
    x2 => plt2%get_x_axis()
    y2 => plt2%get_y_axis()
    lgnd => plt1%get_legend()

    call x1%set_title("Frequency [Hz]")
    call y1%set_title("Amplitude [dB]")

    call x2%set_title("Frequency [Hz]")
    call y2%set_title("Phase [deg]")

    call lgnd%set_is_visible(.true.)
    call lgnd%set_horizontal_position(LEGEND_CENTER)

    call pd1%define_data(freq, amp)
    call pd1%set_name("Actual")
    call pd1%set_line_width(2.0)
    call plt1%push(pd1)

    call pd2%define_data(freq, mdlamp)
    call pd2%set_name("Fitted")
    call pd2%set_line_width(2.0)
    call pd2%set_line_style(LINE_DASHED)
    call plt1%push(pd2)

    call pd1%define_data(freq, phase)
    call plt2%push(pd1)

    call pd2%define_data(freq, mdlphase)
    call plt2%push(pd2)

    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end program