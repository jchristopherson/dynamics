! This example analyzes a planar four-bar linkage.  The mechanism contains a
! single closed kinematic loop, so its forward kinematics require the solution
! of the loop-closure constraints rather than a simple accumulation of link
! transformations.

program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use linalg, only : identity
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    integer(int32), parameter :: npts = 100

    ! Model Properties
    real(real64), parameter :: crank = 1.0d0
    real(real64), parameter :: coupler = 3.5d0
    real(real64), parameter :: rocker = 3.0d0
    real(real64), parameter :: ground = 4.0d0

    ! Local Variables
    integer(int32) :: i
    type(link_container) :: links(4)
    type(joint) :: joints(4)
    type(planar_linkage) :: linkage
    real(real64) :: theta(npts), p(npts,3)
    real(real64), allocatable, dimension(:) :: plt_config
    real(real64), allocatable, dimension(:,:) :: jac

    ! Define the links.  Each link carries a joint frame at each end, with the
    ! body frame located at the first of the two.
    allocate(links(1)%item, source = planar_link(ground))
    allocate(links(2)%item, source = planar_link(crank))
    allocate(links(3)%item, source = planar_link(coupler))
    allocate(links(4)%item, source = planar_link(rocker))

    ! Connect the links.  The crank is the driven member.
    joints(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, actuated = .true.)
    joints(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
    joints(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
    joints(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)

    ! Build the mechanism.  The end-effector is the distal end of the coupler.
    linkage = planar_linkage(links, joints, base = 1, effector = 3, &
        tool = translation(coupler))

    ! Display the mobility of the mechanism
    print "(A,I0)", "Number of independent loops: ", linkage%get_loop_count()
    print "(A,I0)", "Number of joint variables: ", linkage%get_variable_count()
    print "(A,I0)", "Number of constraint equations: ", &
        linkage%get_constraint_count()
    print "(A,I0)", "Degrees of freedom: ", linkage%get_degrees_of_freedom()

    ! Establish a starting configuration.  A closed-loop mechanism admits more
    ! than one assembly mode, so the starting estimate selects the branch.
    call linkage%set_configuration([0.0d0, 0.5d0 * pi, -0.5d0 * pi, 0.0d0])

    ! Sweep the crank and trace the coupler point
    theta = linspace(0.0d0, 2.0d0 * pi, npts)
    do i = 1, npts
        p(i,:) = linkage%end_effector_pose([theta(i)])
    end do

    ! Examine the sensitivity of the coupler point to the crank at mid-stroke
    jac = linkage%jacobian([0.5d0 * pi])
    print "(A)", new_line('a') // "Jacobian at a crank angle of 90 degrees:"
    do i = 1, size(jac, 1)
        print *, jac(i,:)
    end do

    ! Plot the coupler curve along with the linkage. The complete set of joint 
    ! variables, not just the actuated variable, is required in order to locate 
    ! each link.
    plt_config = linkage%solve_configuration([0.25d0 * pi])
    call plot_coupler_path(p, linkage, plt_config)

contains
    ! Constructs a planar link of the requested length carrying a joint frame at
    ! each end.
    function planar_link(length) result(rst)
        real(real64), intent(in) :: length
        type(multi_joint_link) :: rst
        real(real64) :: frames(4, 4, 2)
        frames(:,:,1) = translation(0.0d0)
        frames(:,:,2) = translation(length)
        rst = multi_joint_link(frames)
    end function

    ! A 4-by-4 translation along the x-axis.
    function translation(a) result(rst)
        real(real64), intent(in) :: a
        real(real64) :: rst(4, 4)
        rst = identity(4)
        rst(1,4) = a
    end function

    ! Coupler Path Plot
    subroutine plot_coupler_path(pts, mech, q)
        real(real64), intent(in), dimension(:,:) :: pts
        class(kinematic_mechanism), intent(in), target :: mech
        real(real64), intent(in), dimension(:) :: q

        type(plot_2d) :: plt
        type(legend), pointer :: lgnd

        call plt%initialize()
        call plt%set_x_axis_title("x")
        call plt%set_y_axis_title("y")
        lgnd => plt%get_legend()
        call lgnd%set_is_visible(.true.)
        call lgnd%set_horizontal_position(LEGEND_LEFT)
        call lgnd%set_vertical_position(LEGEND_BOTTOM)
        call lgnd%set_draw_border(.false.)
        call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
        call lgnd%set_draw_inside_axes(.false.)

        ! Coupler Curve
        call plt%push(pts(:,1), pts(:,2), name = "Coupler Curve")

        ! The linkage in its requested configuration
        call draw_linkage(plt, mech, q)

        call plt%draw()
    end subroutine

    ! Draws each link as a polyline passing through the joint frames the link
    ! carries.  The location of each frame is found by transforming the frame,
    ! which is expressed in the link's body coordinate frame, into the frame of
    ! the base link.
    subroutine draw_linkage(plt, mech, q)
        type(plot_2d), intent(inout) :: plt
        class(kinematic_mechanism), intent(in), target :: mech
        real(real64), intent(in), dimension(:) :: q

        integer(int32) :: j, k, n
        class(link), pointer :: lnk
        type(plot_data_2d), allocatable, dimension(:) :: pd
        real(real64) :: T(4, 4), Pk(4, 4)
        real(real64), allocatable, dimension(:) :: x, y
        character(len = :), allocatable :: name

        allocate(pd(mech%get_link_count()))
        do j = 1, mech%get_link_count()
            ! Get the link and compute its location and orientation
            lnk => mech%get_link(j)
            n = lnk%get_joint_count()
            T = mech%body_transform(j, q)
            allocate(x(n), y(n))
            do k = 1, n
                Pk = matmul(T, lnk%get_joint_frame(k))
                x(k) = Pk(1,4)
                y(k) = Pk(2,4)
            end do

            ! Define a name for the link
            select case (j)
            case (1)
                name = "Ground Link"
            case (2)
                name = "Crank"
            case (3)
                name = "Coupler"
            case (4)
                name = "Rocker"
            end select

            ! Plot the link
            call pd(j)%define_data(x, y)
            call pd(j)%set_line_width(3.0)
            if (j == 1) call pd(j)%set_line_style(LINE_DASHED)
            call pd(j)%set_draw_markers(.true.)
            call pd(j)%set_marker_style(MARKER_EMPTY_CIRCLE)
            call pd(j)%set_name(name)
            call plt%push(pd(j))

            deallocate(x, y)
        end do
    end subroutine
end program
