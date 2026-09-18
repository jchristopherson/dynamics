! This example simulates a motor-driven planar four-bar linkage. The crank,
! coupler, and rocker have distributed mass properties, while the ground link
! remains fixed. A sinusoidal angular displacement is prescribed at the crank
! input, and the linkage_dynamic_model computes the required motor torque while
! enforcing every revolute joint and loop constraint.

program example
    use iso_fortran_env, only : int32, real64
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
    integer(int32) :: i
    real(real64), dimension(4) :: q
    real(real64), dimension(2) :: closure_point
    real(real64), dimension(3,3) :: rotation
    real(real64), allocatable, dimension(:) :: time, crank_angle, &
        crank_target, coupler_angle, rocker_angle, motor_torque, residual
    real(real64), allocatable, dimension(:,:) :: constraint_multipliers
    type(link_container), dimension(4) :: links
    type(joint), dimension(4) :: joints
    type(planar_linkage) :: mechanism
    type(linkage_dynamic_model) :: dynamic_model
    type(variational_integrator) :: integrator
    type(variational_state), allocatable, dimension(:) :: solution
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2, plt3, plt4

    ! Define the fixed ground and the three moving massive links. Each moving
    ! link uses a rectangular-prism inertia tensor about its center of mass.
    allocate(links(1)%item, source = massive_planar_link(ground_length, 1.0d0))
    allocate(links(2)%item, source = massive_planar_link(crank_length, &
        crank_mass))
    allocate(links(3)%item, source = massive_planar_link(coupler_length, &
        coupler_mass))
    allocate(links(4)%item, source = massive_planar_link(rocker_length, &
        rocker_mass))

    ! Connect the links with revolute joints. The first joint identifies the
    ! crank as the actuated kinematic input.
    joints(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, actuated = .true.)
    joints(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
    joints(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
    joints(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)
    mechanism = planar_linkage(links, joints, base = 1, effector = 3)

    ! Select an assembly mode and construct the maximal-coordinate dynamic
    ! model from the compatible joint configuration.
    call four_bar_configuration(initial_crank_angle, q, closure_point)
    call mechanism%set_configuration(q)
    dynamic_model = linkage_dynamic_model(mechanism, q)

        ! The base link is omitted from the moving-body array, so body 1 is the
        ! crank. Prescribing its absolute angle adds one rheonomic constraint. The
        ! final constraint multiplier is the motor torque required to follow it.
    integrator%settings%linear_solver = VI_DENSE_SOLVER
    solution = dynamic_model%solve(integrator, dt, ntime, &
		gravity = [0.0d0, -9.80665d0, 0.0d0], &
                prescribed_body = 1, prescribed_motion = crank_motion, &
                multipliers = constraint_multipliers)

    ! Recover the absolute orientation of each moving link.
    allocate(time(ntime), crank_angle(ntime), crank_target(ntime), &
        coupler_angle(ntime), rocker_angle(ntime), motor_torque(ntime))
    do i = 1, ntime
        time(i) = solution(i)%time
        crank_target(i) = crank_motion(time(i)) * 180.0d0 / pi
        rotation = solution(i)%orientation(1)%to_matrix()
        crank_angle(i) = atan2(rotation(2,1), rotation(1,1)) * 180.0d0 / pi
        rotation = solution(i)%orientation(2)%to_matrix()
        coupler_angle(i) = atan2(rotation(2,1), rotation(1,1)) * 180.0d0 / pi
        rotation = solution(i)%orientation(3)%to_matrix()
        rocker_angle(i) = atan2(rotation(2,1), rotation(1,1)) * 180.0d0 / pi
    end do
    motor_torque(2:ntime) = constraint_multipliers( &
        size(constraint_multipliers,1),:)
    motor_torque(1) = motor_torque(2)

    residual = dynamic_model%constraint_residual(solution(ntime))
    print "(A,F7.3,A)", "Simulated ", time(ntime), " seconds."
    print "(A,F9.3,A)", "Final crank angle: ", crank_angle(ntime), " deg"
    print "(A,F9.3,A)", "Peak motor torque: ", &
        maxval(abs(motor_torque)), " N m"
    print "(A,ES10.3)", "Maximum final constraint residual: ", &
        maxval(abs(residual))

    ! Plot all moving-link angles to illustrate the coupled four-bar response.
    call plt%initialize(4, 1, width = 1400, height = 900)
    call plt%set_font_size(12)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()

    call plt1%set_title("Crank")
    call plt1%set_x_axis_title("Time [s]")
    call plt1%set_y_axis_title("Absolute Link Angle [deg]")
    call plt1%push(time, crank_angle, name = "Crank")
    call plt%set(1, 1, plt1)

    call plt2%set_title("Coupler")
    call plt2%set_x_axis_title("Time [s]")
    call plt2%set_y_axis_title("Absolute Link Angle [deg]")
    call plt2%push(time, coupler_angle, name = "Coupler")
    call plt%set(2, 1, plt2)

    call plt3%set_title("Rocker")
    call plt3%set_x_axis_title("Time [s]")
    call plt3%set_y_axis_title("Absolute Link Angle [deg]")
    call plt3%push(time, rocker_angle, name = "Rocker")
    call plt%set(3, 1, plt3)

    call plt4%set_title("Required Motor Torque")
    call plt4%set_x_axis_title("Time [s]")
    call plt4%set_y_axis_title("Torque [N m]")
    call plt4%push(time, motor_torque, name = "Motor Torque")
    call plt%set(4, 1, plt4)

    call plt%draw()

contains
! ------------------------------------------------------------------------------
pure function crank_motion(t) result(rst)
    !! Prescribes one cycle of sinusoidal crank motion centered at 40 degrees
    !! with an amplitude of 10 degrees. The cosine phase starts the mechanism
    !! at its lower displacement limit with zero commanded velocity.
    real(real64), intent(in) :: t
        !! The simulation time.
    real(real64) :: rst
        !! The prescribed absolute crank angle, in radians.

    rst = motion_center - motion_amplitude * &
        cos(2.0d0 * pi * motion_frequency * t)
end function

! ------------------------------------------------------------------------------
function massive_planar_link(length, mass) result(rst)
    !! Constructs a massive rectangular link with joint frames at both ends.
    real(real64), intent(in) :: length
        !! The distance between the two joint centers.
    real(real64), intent(in) :: mass
        !! The link mass.
    type(multi_joint_link) :: rst
        !! The resulting two-joint link.

    real(real64), dimension(4,4,2) :: frames
    real(real64), dimension(3,3) :: inertia

    frames(:,:,1) = translation(0.0d0)
    frames(:,:,2) = translation(length)
    inertia = 0.0d0
    inertia(1,1) = mass * link_width**2 / 6.0d0
    inertia(2,2) = mass * (length**2 + link_width**2) / 12.0d0
    inertia(3,3) = inertia(2,2)
    rst = multi_joint_link(frames, mass = mass, inertia = inertia, &
        cg = [0.5d0 * length, 0.0d0, 0.0d0])
end function

! ------------------------------------------------------------------------------
pure function translation(distance) result(rst)
    !! Constructs a homogeneous translation along the x-axis.
    real(real64), intent(in) :: distance
        !! The translation distance.
    real(real64), dimension(4,4) :: rst
        !! The homogeneous transformation matrix.
    integer(int32) :: i

    rst = 0.0d0
    do i = 1, 4
        rst(i,i) = 1.0d0
    end do
    rst(1,4) = distance
end function

! ------------------------------------------------------------------------------
subroutine four_bar_configuration(theta, q, point)
    !! Computes one compatible assembly of the four-bar linkage for a supplied
    !! crank angle by intersecting the coupler and rocker circles.
    real(real64), intent(in) :: theta
        !! The crank angle.
    real(real64), intent(out), dimension(4) :: q
        !! The four relative revolute-joint angles.
    real(real64), intent(out), dimension(2) :: point
        !! The coupler-rocker joint location.

    real(real64), dimension(2) :: crank_tip, fixed_pivot, direction, normal
    real(real64) :: distance, along, height, coupler_orientation, &
        rocker_orientation

    crank_tip = crank_length * [cos(theta), sin(theta)]
    fixed_pivot = [ground_length, 0.0d0]
    distance = norm2(fixed_pivot - crank_tip)
    direction = (fixed_pivot - crank_tip) / distance
    normal = [-direction(2), direction(1)]
    along = 0.5d0 * (distance**2 + coupler_length**2 - &
        rocker_length**2) / distance
    height = sqrt(coupler_length**2 - along**2)
    point = crank_tip + along * direction + height * normal
    coupler_orientation = atan2(point(2) - crank_tip(2), &
        point(1) - crank_tip(1))
    rocker_orientation = atan2(point(2), point(1) - ground_length)
    q = [theta, coupler_orientation - theta, &
        rocker_orientation - coupler_orientation, -rocker_orientation]
end subroutine

end program