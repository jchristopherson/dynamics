! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
! This example analyzes a planar four-bar linkage.  The mechanism contains a
! single closed kinematic loop, so its forward kinematics require the solution
! of the loop-closure constraints rather than a simple accumulation of link
! transformations.
program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    implicit none

    ! Simulation Parameters
    integer(int32), parameter :: ntime = 1001
    real(real64), parameter :: dt = 5.0d-4
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: motion_frequency = 5.0d0
    real(real64), parameter :: motion_center = 40.0d0 * pi / 180.0d0
    real(real64), parameter :: motion_amplitude = 20.0d0 * pi / 180.0d0
    real(real64), parameter :: torsional_spring_rate = 5.0d3
    real(real64), parameter :: torsional_damping_rate = 1.5d2
    real(real64), parameter :: gc = 9.81d0

    ! Four-Bar Geometry and Mass Properties
    real(real64), parameter :: crank_length = 1.0d0
    real(real64), parameter :: coupler_length = 3.5d0
    real(real64), parameter :: rocker_length = 3.0d0
    real(real64), parameter :: ground_length = 4.0d0
    real(real64), parameter :: crank_mass = 1.0d0
    real(real64), parameter :: coupler_mass = 1.5d0
    real(real64), parameter :: rocker_mass = 1.2d0
    real(real64), parameter :: link_width = 5.0d-2
    real(real64), parameter :: initial_crank_angle = &
        motion_center - motion_amplitude

    ! Local Variables
    type(link_container) :: links(4)
    type(joint) :: joints(4)
    type(planar_linkage) :: mechanism
    type(linkage_dynamic_model) :: model
    procedure(linkage_prescribed_motion), pointer :: motor_motion
    type(torsional_spring) :: spring
    type(torsional_damper) :: damper
    real(real64) :: g(3), R(3,3)
    type(variational_integrator) :: integrator
    type(variational_state), allocatable, dimension(:) :: solution
    real(real64), allocatable, dimension(:,:) :: constraint_multipliers
    type(joint_reaction) :: joint_reactions(4)
    integer(int32) :: i, nmult
    real(real64), allocatable, dimension(:) :: time, torque, crank_angle, &
        rocker_angle, Fx1, Fy1, Fx2, Fy2, Fx3, Fy3, Fx4, Fy4

    ! Plot Variables
    type(multiplot) :: plt, fplt
    type(plot_2d) :: plt1, plt2, plt3, fplt1, fplt2, fplt3, fplt4
    class(legend), pointer :: lgnd

! ******************************************************************************
! BUILD THE LINKAGE
! ------------------------------------------------------------------------------
    ! The list of links
    allocate(links(1)%item, source = create_link(ground_length, link_width, 1.0d0))
    allocate(links(2)%item, source = create_link(crank_length, link_width, crank_mass))
    allocate(links(3)%item, source = create_link(coupler_length, link_width, coupler_mass))
    allocate(links(4)%item, source = create_link(rocker_length, link_width, rocker_mass))

    ! Connect the links via revolute joints
    joints(1) = joint( &
        REVOLUTE_JOINT, &       ! joint type
        1, &                    ! parent link index (ground link)
        2, &                    ! child link index (crank link)
        1, &                    ! parent coordinate frame index (interfaces to joint #1 on the ground link)
        1, &                    ! child coordinate frame index  (interfaces to joint #1 on the crank link)
        actuated = .true. &     ! we're going to drive this joint with a motor
    )
    joints(2) = joint( &
        REVOLUTE_JOINT, &       ! joint type
        2, &                    ! parent link index (crank link)
        3, &                    ! child link index (coupler link)
        2, &                    ! parent coordinate frame index (interfaces to joint #2 on the crank)
        1 &                     ! child coordinate frame index (interfaces to joint #1 on the coupler)
    )
    joints(3) = joint( &
        REVOLUTE_JOINT, &       ! joint type
        3, &                    ! parent link index (coupler link)
        4, &                    ! child link index (rocker link)
        2, &                    ! parent coordinate frame index (interfaces to joint #2 on the coupler)
        2 &                     ! child coordinate frame index (interfaces to joint #2 on the rocker)
    )
    joints(4) = joint( &
        REVOLUTE_JOINT, &       ! joint type
        4, &                    ! parent link index (rocker link)
        1, &                    ! child link index (ground link)
        1, &                    ! parent coordinate frame index (interfaces to joint #1 on the rocker)
        2 &                     ! child coordinate frame index (interfaces to joint #2 on the ground link)
    )

    ! Construct the mechanism
    mechanism = planar_linkage( &
        links, &                ! the list of link objects
        joints, &               ! the list of joint objects
        base = 1, &             ! the index of the base link (ground)
        effector = 3 &          ! the index of the end-effector link (coupler)
    )

! ******************************************************************************
! SET UP THE DYNAMIC ANALYSIS
! ------------------------------------------------------------------------------
    ! As a closed-loop mechanism can have multiple valid configurations, we 
    ! need to establish a reasonably close estimate to the initial geometry
    ! to ensure the solution tracks on the correct configuration.
    call define_initial_configuration(mechanism, initial_crank_angle)

    ! Define the dynamic model
    model = linkage_dynamic_model( &
        mechanism, &                        ! the mechanism
        mechanism%get_configuration() &     ! the initial configuration
    )

    ! Now we can define the necessary routines for the integrator
    motor_motion => crank_motor

    ! Add a spring and damper element
    spring%free_angle = 0.0d0
    spring%joint_index = 4      ! the spring is tied to the rocker-ground connection
    spring%stiffness = torsional_spring_rate
    call model%add_torsional_spring(spring)

    damper%joint_index = 4      ! the damper is tied to the rocker-ground connection
    damper%damping = torsional_damping_rate
    call model%add_torsional_damper(damper)

    ! Define the gravitational vector
    g = [0.0d0, -gc, 0.0d0]     ! gravity in -y direction

    ! Set up the integrator & solve
    integrator%settings%linear_solver = VI_DENSE_SOLVER ! define the solver type
    solution = model%solve( &
        integrator, &                           ! the integrator to utilize
        dt, &                                   ! time step size
        ntime, &                                ! # of time steps
        gravity = g, &                          ! gravitational vector
        prescribed_body = 1, &                  ! the index of the driven body - ground is ignored in this instance, so the driving body is the crank and thus index 1
        prescribed_motion = motor_motion, &     ! motor motion routine
        multipliers = constraint_multipliers &  ! Lagrange multiplier values - use to get motor torque
    )

! ******************************************************************************
! EXTRACT THE SOLUTION PARAMETERS OF INTEREST
! ------------------------------------------------------------------------------
    allocate( &
        time(ntime), &
        torque(ntime), &
        crank_angle(ntime), &
        rocker_angle(ntime), &
        Fx1(ntime), Fy1(ntime), &
        Fx2(ntime), Fy2(ntime), &
        Fx3(ntime), Fy3(ntime), &
        Fx4(ntime), Fy4(ntime) &
    )
    do i = 1, ntime
        time(i) = solution(i)%time

        ! Crank Angle
        R = solution(i)%orientation(1)%to_matrix()
        crank_angle(i) = atan2(R(2,1), R(1,1)) * 1.8d2 / pi

        ! Rocker Pivot Angle
        R = solution(i)%orientation(3)%to_matrix()
        rocker_angle(i) = atan2(R(2,1), R(1,1)) * 1.8d2 / pi

        ! Joint Reactions
        joint_reactions = model%get_joint_reactions( &
            solution(i), &
            constraint_multipliers(:,i) &
        )
        
        Fx1(i) = joint_reactions(1)%force(1)
        Fy1(i) = joint_reactions(1)%force(2)

        Fx2(i) = joint_reactions(2)%force(1)
        Fy2(i) = joint_reactions(2)%force(2)

        Fx3(i) = joint_reactions(3)%force(1)
        Fy3(i) = joint_reactions(3)%force(2)

        Fx4(i) = joint_reactions(4)%force(1)
        Fy4(i) = joint_reactions(4)%force(2)
    end do

    ! Get the motor torque
    nmult = size(constraint_multipliers, 1)
    torque = constraint_multipliers(nmult,:)    ! motor torque

! ******************************************************************************
! PLOT THE SOLUTION PARAMETERS
! ------------------------------------------------------------------------------
    call plt%initialize(3, 1, width = 1200, height = 800)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()

    call plt1%push(time, crank_angle)
    call plt1%set_x_axis_title("t [sec]")
    call plt1%set_y_axis_title("{/Symbol q}_1 [deg]")
    call plt%set(1, 1, plt1)

    call plt2%push(time, rocker_angle)
    call plt2%set_x_axis_title("t [sec]")
    call plt2%set_y_axis_title("{/Symbol q}_4 [deg]")
    call plt%set(2, 1, plt2)

    call plt3%push(time, torque)
    call plt3%set_x_axis_title("t [sec]")
    call plt3%set_y_axis_title("T [N*m]")
    call plt%set(3, 1, plt3)

    call plt%draw()

! ------------------------------------------------------------------------------
    ! Plot the joint loads
    call fplt%initialize(2, 2, width = 1400, height = 800)
    call fplt1%initialize()
    call fplt2%initialize()
    call fplt3%initialize()
    call fplt4%initialize()
    lgnd => fplt1%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    call lgnd%set_is_opaque(.false.)

    call fplt1%push(time, Fx1, name = "F_x")
    call fplt1%push(time, Fy1, name = "F_y")
    call fplt1%set_x_axis_title("t [sec]")
    call fplt1%set_y_axis_title("F_1 [N]")
    call fplt%set(1, 1, fplt1)

    call fplt2%push(time, Fx2, name = "F_x")
    call fplt2%push(time, Fy2, name = "F_y")
    call fplt2%set_x_axis_title("t [sec]")
    call fplt2%set_y_axis_title("F_2 [N]")
    call fplt%set(2, 1, fplt2)

    call fplt3%push(time, Fx3, name = "F_x")
    call fplt3%push(time, Fy3, name = "F_y")
    call fplt3%set_x_axis_title("t [sec]")
    call fplt3%set_y_axis_title("F_3 [N]")
    call fplt%set(1, 2, fplt3)

    call fplt4%push(time, Fx4, name = "F_x")
    call fplt4%push(time, Fy4, name = "F_y")
    call fplt4%set_x_axis_title("t [sec]")
    call fplt4%set_y_axis_title("F_4 [N]")
    call fplt%set(2, 2, fplt4)

    call fplt%draw()
contains
! ------------------------------------------------------------------------------
    function crank_motor(t, args) result(rst)
        !! Defines the position of the crank motor at time t.
        real(real64), intent(in) :: t
            !! The simulation time.
        class(*), intent(inout), optional :: args
            !! Optional user-defined communication argument.
        real(real64) :: rst
            !! The motor position, in radians.

        rst = motion_center - motion_amplitude * cos(2.0d0 * pi * motion_frequency * t)
    end function

! ------------------------------------------------------------------------------
    function create_link(length_, width_, mass_) result(rst)
        !! Constructs a new 2-joint link with mass properties.
        real(real64), intent(in) :: length_
            !! The length of the link.
        real(real64), intent(in) :: width_
            !! The width of the link.
        real(real64), intent(in) :: mass_
            !! The mass of the link.
        type(multi_joint_link) :: rst
            !! The link.

        ! Variables
        real(real64) :: frames(4, 4, 2) ! 2x 4-by-4 homogeneous transformation matrices
        real(real64) :: inertia(3, 3)   ! 3-by-3 inertia matrix

        ! Define the joint locations
        frames(:,:,1) = translate(0.0d0, 0.0d0, 0.0d0)
        frames(:,:,2) = translate(length_, 0.0d0, 0.0d0)

        ! Define the inertia information
        inertia = 0.0d0
        inertia(1,1) = mass_ * width_**2 / 6.0d0
        inertia(2,2) = mass_ * (length_**2 + width_**2) / 12.0d0
        inertia(3,3) = inertia(2,2)
        rst = multi_joint_link(frames, mass = mass_, inertia = inertia, & 
            cg = [0.5d0 * length_, 0.0d0, 0.0d0])
    end function

! ------------------------------------------------------------------------------
    subroutine define_initial_configuration(mech_, angle_)
        !! Defines a valid configuration estimate for the mechanism.
        class(planar_linkage), intent(inout) :: mech_
            !! The linkage.
        real(real64), intent(in) :: angle_
            !! The initial crank angle, in radians.

        ! Variables
        real(real64) :: cl, cx, cy, rx, ry, rl, dist, dx, dy, nx, ny, cpl, &
            along, height, px, py, cpa, ra, gl, q(4)

        ! Determine the location of joint 2 on the crank
        cl = crank_length
        cx = cl * cos(angle_)   ! x location of joint 2 of the crank
        cy = cl * sin(angle_)   ! y location of joint 2 of the crank

        ! Determine the position of the rocker pivot
        rx = ground_length
        ry = 0.0d0

        ! Determine the ground length
        gl = ground_length

        ! Get the rocker length
        rl = rocker_length

        ! Determine the coupler link length
        cpl = coupler_length
        
        ! Locate the connection point between the rocker and coupler
        dx = rx - cx
        dy = ry - cy
        dist = sqrt(dx**2 + dy**2)
        dx = dx / dist
        dy = dy / dist
        nx = -dy
        ny = dx
        along = 0.5d0 * (dist**2 + cpl**2 - rl**2) / dist
        height = sqrt(cpl**2 - along**2)
        px = cx + along * dx + height * nx  ! x-coordinate of rocker end point
        py = cy + along * dy + height * ny  ! y-coordinate of rocker end point

        ! Determine the orientation angles
        cpa = atan2(py - cy, px - cx)
        ra = atan2(py, px - gl)

        ! Define the initial joint variables
        q = [angle_, cpa - angle_, ra - cpa, -ra]
        call mech_%set_configuration(q)
    end subroutine

! ------------------------------------------------------------------------------
end program