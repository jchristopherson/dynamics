program example
    use dynamics
    use fplot_core
    use iso_fortran_env
    implicit none

    ! Parameters
    real(real64), parameter :: r1(3) = [1.0d0, 1.0d0, 0.0d0]
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: x_angle = pi / 4.0d0
    real(real64), parameter :: z_angle = pi / 4.0d0

    ! Local Variables
    real(real64) :: Rx(3, 3), Rz(3, 3), Rzx(3, 3), r2a(3), r3a(3), r3b(3)

    ! Plot Variables
    type(plot_3d) :: plt
    class(plot_axis), pointer :: xAxis, yAxis, zAxis
    class(legend), pointer :: lgnd
    type(plot_data_3d) :: d1, d2a, d3a, d3b
    
    ! Compute rotation matrices
    Rx = rotate_x(x_angle)
    Rz = rotate_z(z_angle)
    Rzx = matmul(Rx, Rz)    ! rotate about z, and then about x

    ! Compute the rotation vectors individually
    r2a = matmul(Rz, r1)
    r3a = matmul(Rx, r2a)

    ! Now, as a whole
    r3b = matmul(Rzx, r1)

    ! Show a plot illustrating the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    zAxis => plt%get_z_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("x")
    call yAxis%set_title("y")
    call zAxis%set_title("z")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_LEFT)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_draw_border(.false.)
    call plt%set_axis_equal(.true.)

    call d1%define_data( &
        [0.0d0, r1(1)], &
        [0.0d0, r1(2)], &
        [0.0d0, r1(3)] &
    )
    call d1%set_line_width(2.0)
    call d1%set_name("r_1")
    call plt%push(d1)

    call d2a%define_data( &
        [0.0d0, r2a(1)], &
        [0.0d0, r2a(2)], &
        [0.0d0, r2a(3)] &
    )
    call d2a%set_line_width(2.0)
    call d2a%set_line_style(LINE_DASHED)
    call d2a%set_name("R_z r_1")
    call plt%push(d2a)

    call d3a%define_data( &
        [0.0d0, r3a(1)], &
        [0.0d0, r3a(2)], &
        [0.0d0, r3a(3)] &
    )
    call d3a%set_line_width(2.0)
    call d3a%set_line_style(LINE_DASH_DOTTED)
    call d3a%set_name("R_x R_z r_1")
    call plt%push(d3a)

    call d3b%define_data( &
        [0.0d0, r3b(1)], &
        [0.0d0, r3b(2)], &
        [0.0d0, r3b(3)] &
    )
    call d3b%set_line_width(2.0)
    call d3b%set_line_style(LINE_DASH_DOT_DOT)
    call d3b%set_name("R_{zx} r_1")
    call plt%push(d3b)

    call plt%draw()
end program