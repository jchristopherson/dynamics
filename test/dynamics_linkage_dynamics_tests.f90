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
    real(real64), allocatable, dimension(:) :: residual
    type(binary_link), dimension(1) :: links
    type(serial_linkage) :: mechanism
    type(linkage_dynamic_model) :: model
    type(variational_integrator) :: integrator
    type(variational_state) :: initial
    type(variational_state), allocatable, dimension(:) :: solution

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
end function

! ------------------------------------------------------------------------------
function test_parallel_linkage_dynamics() result(rst)
    !! Verifies maximal-coordinate conversion of a closed planar four-bar.
    logical :: rst
    real(real64), parameter :: crank = 1.0d0
    real(real64), parameter :: coupler = 3.5d0
    real(real64), parameter :: rocker = 3.0d0
    real(real64), parameter :: ground = 4.0d0
    real(real64), dimension(4) :: q
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