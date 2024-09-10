module excitation
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
    use excitation
    implicit none

    ! Parameters
    integer(int32), parameter :: nfreq = 1000
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: fmin = 2.0d0 * pi * 10.0d0
    real(real64), parameter :: fmax = 2.0d0 * pi * 1.0d3
    real(real64), parameter :: alpha = 1.0d-3
    real(real64), parameter :: beta = 2.0d-6

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
    real(real64), allocatable, dimension(:) :: freq
    real(real64), allocatable, dimension(:,:) :: mag, phase
    procedure(modal_excite), pointer :: fcn

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd1, pd2, pd3
    class(plot_axis), pointer :: xAxis, yAxis
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

    ! Extract the magnitude and phase information.
    freq = rsp%frequency / (2.0d0 * pi)
    mag = abs(rsp%responses)
    phase = (1.8d2 / pi) * atan2(aimag(rsp%responses), real(rsp%responses))

! ------------------------------------------------------------------------------
! PLOTTING CODE ONLY
! ------------------------------------------------------------------------------
    ! Plot the frequency response functions
    call plt%initialize(2, 1)
    call plt1%initialize()
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    lgnd => plt1%get_legend()
    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("|X|")
    call yAxis%set_is_log_scaled(.true.)
    call yAxis%set_use_default_tic_label_format(.false.)
    call yAxis%set_tic_label_format("%0.0e")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    call lgnd%set_draw_border(.false.)

    call pd1%define_data(freq, mag(:,1))
    call pd1%set_line_width(2.0)
    call pd1%set_name("X_1")
    call plt1%push(pd1)

    call pd2%define_data(freq, mag(:,2))
    call pd2%set_line_width(2.0)
    call pd2%set_line_style(LINE_DASHED)
    call pd2%set_name("X_2")
    call plt1%push(pd2)

    call pd3%define_data(freq, mag(:,3))
    call pd3%set_line_width(3.0)
    call pd3%set_line_style(LINE_DASH_DOTTED)
    call pd3%set_name("X_3")
    call plt1%push(pd3)

    call plt2%initialize()
    xAxis => plt2%get_x_axis()
    yAxis => plt2%get_y_axis()
    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("{/Symbol f} [deg]")
    call yAxis%set_use_default_tic_label_format(.false.)
    call yAxis%set_tic_label_format("%4.0f")

    call pd1%define_data(freq, phase(:,1))
    call plt2%push(pd1)

    call pd2%define_data(freq, phase(:,2))
    call plt2%push(pd2)

    call pd3%define_data(freq, phase(:,3))
    call plt2%push(pd3)
    
    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end program