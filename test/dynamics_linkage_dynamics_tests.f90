module dynamics_linkage_dynamics_tests
    use iso_fortran_env, only : int32, real64
    use fortran_test_helper
    use dynamics
    implicit none

contains
! ------------------------------------------------------------------------------
function test_serial_linkage_dynamics() result(rst)
    !! Verifies serial-link conversion, gravity response, and constraint
    !! preservation for a massive revolute link.
    logical :: rst
    real(real64), parameter :: angle = 0.4d0
    real(real64), dimension(3,3) :: inertia
    real(real64), dimension(3,1) :: body_force, body_torque
    real(real64), dimension(5) :: test_multipliers
    real(real64), allocatable, dimension(:) :: residual
    type(binary_link), dimension(1) :: links
    type(serial_linkage) :: mechanism
    type(linkage_dynamic_model) :: model
    type(variational_integrator) :: integrator
    type(variational_state) :: initial
    type(variational_state), allocatable, dimension(:) :: solution
    type(joint_reaction), allocatable, dimension(:) :: reactions

    rst = .true.
    inertia = 0.0d0
    inertia(1,1) = 1.0d0 / 12.0d0
    inertia(2,2) = 1.0d-3
    inertia(3,3) = 1.0d0 / 12.0d0
    links(1) = binary_link(length = 1.0d0, mass = 1.0d0, &
        inertia = inertia, cg = [-0.5d0, 0.0d0, 0.0d0])
    mechanism = serial_linkage(links)
    model = linkage_dynamic_model(mechanism, [angle])
    initial = model%get_initial_state()

    if (.not.assert(model%get_body_count(), 1)) rst = .false.
    if (.not.assert(model%get_constraint_count(), 5)) rst = .false.
    residual = model%constraint_residual(initial)
    if (maxval(abs(residual)) > 1.0d-10) then
        rst = .false.
        print "(A)", "TEST FAILED: test_serial_linkage_dynamics - initial constraints"
    end if

    test_multipliers = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
    reactions = model%get_joint_reactions(initial, test_multipliers)
    if (.not.assert(reactions(1)%force, [1.0d0, 2.0d0, 3.0d0])) then
        rst = .false.
        print "(A)", "TEST FAILED: test_serial_linkage_dynamics - joint force"
    end if
    if (.not.assert(reactions(1)%moment, [5.0d0, -4.0d0, 0.0d0])) then
        rst = .false.
        print "(A)", "TEST FAILED: test_serial_linkage_dynamics - joint moment"
    end if

    solution = model%solve(integrator, 1.0d-3, 3, &
        gravity = [0.0d0, -9.81d0, 0.0d0])
    residual = model%constraint_residual(solution(3))
    if (maxval(abs(residual)) > 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_serial_linkage_dynamics - final constraints"
    end if
    if (abs(solution(3)%angular_velocity(3,1)) <= 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_serial_linkage_dynamics - gravity response"
    end if

    body_force(:,1) = [1.0d0, 2.0d0, 3.0d0]
    body_torque(:,1) = [0.0d0, 0.0d0, 1.0d0]
    solution = model%solve(integrator, 1.0d-3, 2, &
        initial_state = initial, body_force = body_force, &
        body_torque = body_torque)
    if (abs(solution(2)%angular_velocity(3,1)) <= 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_serial_linkage_dynamics - body torque"
    end if
end function

! ------------------------------------------------------------------------------
function test_spatial_joint_dynamics() result(rst)
    !! Exercises compatible configurations and reaction mappings for the
    !! spatial prismatic, cylindrical, universal, and spherical joints.
    logical :: rst
    integer(int32), parameter :: ncase = 4
    integer(int32), dimension(ncase) :: joint_types
    integer(int32) :: i, j, nconstraint
    real(real64), allocatable, dimension(:) :: q, residual, multipliers
    type(parallel_linkage) :: mechanism
    type(linkage_dynamic_model) :: model
    type(variational_state) :: state
    type(joint_reaction), allocatable, dimension(:) :: reactions

    rst = .true.
    joint_types = [PRISMATIC_JOINT, CYLINDRICAL_JOINT, &
        UNIVERSAL_JOINT, SPHERICAL_JOINT]
    do i = 1, ncase
        q = spatial_joint_configuration(joint_types(i))
        mechanism = build_spatial_joint_mechanism(joint_types(i))
        model = linkage_dynamic_model(mechanism, q)
        state = model%get_initial_state()
        residual = model%constraint_residual(state)
        if (maxval(abs(residual)) > 1.0d-10) then
            rst = .false.
            print "(A,I0)", &
                "TEST FAILED: test_spatial_joint_dynamics - residual ", i
        end if

        nconstraint = model%get_constraint_count()
        allocate(multipliers(nconstraint))
        multipliers = [(real(j,real64), j=1,nconstraint)]
        reactions = model%get_joint_reactions(state, multipliers)
        select case (joint_types(i))
        case (PRISMATIC_JOINT)
            call check_reaction(reactions(1), [1.0d0, 2.0d0, 0.0d0], &
                [3.0d0, 4.0d0, 5.0d0], rst, i)
        case (CYLINDRICAL_JOINT)
            call check_reaction(reactions(1), [1.0d0, 2.0d0, 0.0d0], &
                [4.0d0, -3.0d0, 0.0d0], rst, i)
        case (UNIVERSAL_JOINT)
            call check_reaction(reactions(1), [1.0d0, 2.0d0, 3.0d0], &
                [0.0d0, 0.0d0, -4.0d0], rst, i)
        case (SPHERICAL_JOINT)
            call check_reaction(reactions(1), [1.0d0, 2.0d0, 3.0d0], &
                [0.0d0, 0.0d0, 0.0d0], rst, i)
        end select
        deallocate(multipliers)
    end do
end function

! ------------------------------------------------------------------------------
function test_fixed_and_planar_prismatic_joints() result(rst)
    !! Covers fixed spatial constraints and fixed/prismatic planar constraints,
    !! including their reaction wrench mappings.
    logical :: rst
    integer(int32) :: i
    real(real64), allocatable, dimension(:) :: residual, multipliers
    real(real64), dimension(4,4,1) :: single_frame
    real(real64), dimension(4,4,2) :: double_frame
    type(link_container), dimension(3) :: links
    type(joint), dimension(2) :: joints
    type(parallel_linkage) :: spatial
    type(planar_linkage) :: planar
    type(linkage_dynamic_model) :: model
    type(variational_state) :: state
    type(joint_reaction), allocatable, dimension(:) :: reactions

    rst = .true.

    ! A fixed body followed by a spherical joint exercises the six-component
    ! fixed spatial residual and reaction mapping.
    single_frame(:,:,1) = translate_x(0.0d0)
    allocate(links(1)%item, source = multi_joint_link(single_frame))
    allocate(links(2)%item, source = multi_joint_link(single_frame))
    allocate(links(3)%item, source = multi_joint_link(single_frame))
    joints(1) = joint(FIXED_JOINT, 1, 2)
    joints(2) = joint(SPHERICAL_JOINT, 2, 3)
    spatial = parallel_linkage(links, joints, base = 1)
    model = linkage_dynamic_model(spatial, [0.2d0, -0.1d0, 0.15d0])
    state = model%get_initial_state()
    residual = model%constraint_residual(state)
    if (maxval(abs(residual)) > 1.0d-10) rst = .false.
    multipliers = [(real(i,real64), i=1,model%get_constraint_count())]
    reactions = model%get_joint_reactions(state, multipliers)
    call check_reaction(reactions(1), [1.0d0, 2.0d0, 3.0d0], &
        [4.0d0, 5.0d0, 6.0d0], rst, 5)

    ! Rotate the prismatic joint frame so its z-axis lies along world x. The
    ! moving body may translate along x but is constrained in y and rotation.
    deallocate(links(1)%item, links(2)%item, links(3)%item)
    single_frame(:,:,1) = translate_x(0.0d0)
    double_frame(:,:,1) = translate_x(0.0d0)
    double_frame(:,:,2) = rotate_frame_y(0.5d0 * acos(-1.0d0))
    allocate(links(1)%item, source = multi_joint_link(single_frame))
    allocate(links(2)%item, source = multi_joint_link(double_frame))
    single_frame(:,:,1) = rotate_frame_y(0.5d0 * acos(-1.0d0))
    allocate(links(3)%item, source = multi_joint_link(single_frame))
    joints(1) = joint(FIXED_JOINT, 1, 2, 1, 1)
    joints(2) = joint(PRISMATIC_JOINT, 2, 3, 2, 1)
    planar = planar_linkage(links, joints, base = 1)
    model = linkage_dynamic_model(planar, [0.25d0])
    state = model%get_initial_state()
    residual = model%constraint_residual(state)
    if (maxval(abs(residual)) > 1.0d-10) then
        rst = .false.
        print "(A)", &
            "TEST FAILED: test_fixed_and_planar_prismatic_joints - residual"
    end if
            deallocate(multipliers)
    allocate(multipliers(model%get_constraint_count()))
    multipliers = 0.0d0
    multipliers(7:9) = [1.0d0, 2.0d0, 3.0d0]
    multipliers(10:11) = [4.0d0, 5.0d0]
    reactions = model%get_joint_reactions(state, multipliers)
    call check_reaction(reactions(1), [1.0d0, 2.0d0, 0.0d0], &
        [0.0d0, 0.0d0, 3.0d0], rst, 6)
    call check_reaction(reactions(2), [0.0d0, 4.0d0, 0.0d0], &
        [0.0d0, 0.0d0, 5.0d0], rst, 7)
end function

! ------------------------------------------------------------------------------
pure function rotate_frame_y(angle) result(rst)
    !! Constructs a homogeneous rotation about the y-axis.
    real(real64), intent(in) :: angle
    real(real64), dimension(4,4) :: rst

    rst = translate_x(0.0d0)
    rst(1,1) = cos(angle)
    rst(3,1) = -sin(angle)
    rst(1,3) = sin(angle)
    rst(3,3) = cos(angle)
end function

! ------------------------------------------------------------------------------
function build_spatial_joint_mechanism(joint_type) result(rst)
    !! Builds a ground-to-body mechanism containing one requested joint.
    integer(int32), intent(in) :: joint_type
    type(parallel_linkage) :: rst
    type(link_container), dimension(2) :: links
    type(joint), dimension(1) :: joints
    real(real64), dimension(4,4,1) :: frames

    frames(:,:,1) = translate_x(0.0d0)
    allocate(links(1)%item, source = multi_joint_link(frames))
    allocate(links(2)%item, source = multi_joint_link(frames))
    joints(1) = joint(joint_type, 1, 2, 1, 1, .true.)
    rst = parallel_linkage(links, joints, base = 1)
end function

! ------------------------------------------------------------------------------
function spatial_joint_configuration(joint_type) result(rst)
    !! Supplies a nonzero compatible coordinate vector for a spatial joint.
    integer(int32), intent(in) :: joint_type
    real(real64), allocatable, dimension(:) :: rst

    select case (joint_type)
    case (PRISMATIC_JOINT, REVOLUTE_JOINT)
        rst = [0.2d0]
    case (CYLINDRICAL_JOINT, UNIVERSAL_JOINT)
        rst = [0.2d0, -0.1d0]
    case (SPHERICAL_JOINT)
        rst = [0.2d0, -0.1d0, 0.15d0]
    end select
end function

! ------------------------------------------------------------------------------
subroutine check_reaction(actual, force, moment, rst, test_case)
    !! Compares a computed reaction wrench against expected values.
    type(joint_reaction), intent(in) :: actual
    real(real64), intent(in), dimension(3) :: force, moment
    logical, intent(inout) :: rst
    integer(int32), intent(in) :: test_case

    if (.not.assert(actual%force, force) .or. &
        .not.assert(actual%moment, moment)) then
        rst = .false.
        print "(A,I0)", &
            "TEST FAILED: test_spatial_joint_dynamics - reaction ", test_case
    end if
end subroutine

! ------------------------------------------------------------------------------
function test_parallel_linkage_dynamics() result(rst)
    !! Verifies maximal-coordinate conversion of a closed planar four-bar.
    logical :: rst
    real(real64), parameter :: crank = 1.0d0
    real(real64), parameter :: coupler = 3.5d0
    real(real64), parameter :: rocker = 3.0d0
    real(real64), parameter :: ground = 4.0d0
    real(real64), dimension(4) :: q
    real(real64), dimension(17) :: test_multipliers
    real(real64), dimension(2) :: closure_point
    real(real64), allocatable, dimension(:) :: residual
    real(real64), allocatable, dimension(:,:) :: multipliers
    type(link_container), dimension(4) :: links
    type(joint), dimension(4) :: joints
    type(planar_linkage) :: mechanism
    type(linkage_dynamic_model) :: model
    type(variational_integrator) :: integrator
    type(variational_state) :: state
    type(variational_state), allocatable, dimension(:) :: solution
    type(joint_reaction), allocatable, dimension(:) :: reactions

    rst = .true.
    allocate(links(1)%item, source = planar_dynamic_link(ground, 2.0d0))
    allocate(links(2)%item, source = planar_dynamic_link(crank, 1.0d0))
    allocate(links(3)%item, source = planar_dynamic_link(coupler, 1.5d0))
    allocate(links(4)%item, source = planar_dynamic_link(rocker, 1.2d0))
    joints(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, .true.)
    joints(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
    joints(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
    joints(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)
    mechanism = planar_linkage(links, joints, base = 1)
    call four_bar_configuration(0.7d0, crank, coupler, rocker, ground, q, &
        closure_point)
    model = linkage_dynamic_model(mechanism, q)
    state = model%get_initial_state()

    if (.not.assert(model%get_body_count(), 3)) rst = .false.
    if (.not.assert(model%get_constraint_count(), 17)) rst = .false.
    residual = model%constraint_residual(state)
    if (maxval(abs(residual)) > 1.0d-9) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - constraints"
    end if

    test_multipliers = 0.0d0
    test_multipliers(10:11) = [2.0d0, 3.0d0]
    reactions = model%get_joint_reactions(state, test_multipliers)
    if (.not.assert(model%get_joint_count(), 4)) rst = .false.
    if (.not.assert(reactions(1)%force, [2.0d0, 3.0d0, 0.0d0])) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - joint force"
    end if
    if (.not.assert(reactions(1)%moment, [0.0d0, 0.0d0, 0.0d0])) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - joint moment"
    end if

    solution = model%solve(integrator, 1.0d-4, 2, &
        gravity = [0.0d0, -9.81d0, 0.0d0])
    residual = model%constraint_residual(solution(2))
    if (maxval(abs(residual)) > 1.0d-7) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - dynamic step"
    end if

    solution = model%solve(integrator, 1.0d-4, 2, &
        prescribed_body = 1, prescribed_motion = fixed_crank_motion, &
        multipliers = multipliers)
    if (.not.assert(size(multipliers,1), &
        model%get_constraint_count() + 1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - motor multiplier"
    end if
    if (.not.assert(size(multipliers,2), size(solution))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - multiplier history"
    end if
    reactions = model%get_joint_reactions(solution(2), multipliers(:,1))
    if (.not.assert(size(reactions), model%get_joint_count())) then
        rst = .false.
        print "(A)", "TEST FAILED: test_parallel_linkage_dynamics - reaction count"
    end if
end function

! ------------------------------------------------------------------------------
pure function fixed_crank_motion(t) result(rst)
    !! Supplies the fixed crank angle used to exercise prescribed motion.
    real(real64), intent(in) :: t
    real(real64) :: rst

    rst = 0.7d0
end function

! ------------------------------------------------------------------------------
function planar_dynamic_link(length, mass) result(rst)
    real(real64), intent(in) :: length, mass
    type(multi_joint_link) :: rst
    real(real64), dimension(4,4,2) :: frames
    real(real64), dimension(3,3) :: inertia

    frames(:,:,1) = translate_x(0.0d0)
    frames(:,:,2) = translate_x(length)
    inertia = 0.0d0
    inertia(1,1) = 1.0d-3
    inertia(2,2) = mass * length**2 / 12.0d0
    inertia(3,3) = inertia(2,2)
    rst = multi_joint_link(frames, mass = mass, inertia = inertia, &
        cg = [0.5d0 * length, 0.0d0, 0.0d0])
end function

! ------------------------------------------------------------------------------
pure function translate_x(distance) result(rst)
    real(real64), intent(in) :: distance
    real(real64), dimension(4,4) :: rst
    integer(int32) :: i

    rst = 0.0d0
    do i = 1, 4
        rst(i,i) = 1.0d0
    end do
    rst(1,4) = distance
end function

! ------------------------------------------------------------------------------
subroutine four_bar_configuration(theta, crank, coupler, rocker, ground, q, &
    closure_point)
    real(real64), intent(in) :: theta, crank, coupler, rocker, ground
    real(real64), intent(out), dimension(4) :: q
    real(real64), intent(out), dimension(2) :: closure_point
    real(real64), dimension(2) :: crank_tip, fixed_pivot, direction, normal
    real(real64) :: distance, along, height, coupler_angle, rocker_angle

    crank_tip = [crank * cos(theta), crank * sin(theta)]
    fixed_pivot = [ground, 0.0d0]
    distance = norm2(fixed_pivot - crank_tip)
    direction = (fixed_pivot - crank_tip) / distance
    normal = [-direction(2), direction(1)]
    along = 0.5d0 * (distance**2 + coupler**2 - rocker**2) / distance
    height = sqrt(coupler**2 - along**2)
    closure_point = crank_tip + along * direction + height * normal
    coupler_angle = atan2(closure_point(2) - crank_tip(2), &
        closure_point(1) - crank_tip(1))
    rocker_angle = atan2(closure_point(2), closure_point(1) - ground)
    q = [theta, coupler_angle - theta, rocker_angle - coupler_angle, &
        -rocker_angle]
end subroutine

end module