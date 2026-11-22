program flexible_pendulum_example
    !! Gravity-driven swing and bending of a beam pinned to a fixed stub.
    !! Uses SI units and plots the motion with FPLOT.
    use iso_fortran_env, only : real64
    use dynamics, only : fe_node, degree_of_freedom, finite_element_mesh, euler_bernoulli_beam_element, &
        linear_elastic_isotropic_material, nodal_condition, flexible_joint, REVOLUTE_JOINT, &
        DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ, flexible_transient_settings, solve_flexible_transient, &
        FLEXIBLE_KENNEDY_CARPENTER_5
    use fplot_core, only : multiplot, plot_2d, GNUPLOT_TERMINAL_PNG, CLR_BLUE, CLR_RED
    implicit none

    integer, parameter :: nsamples = 501
    integer, parameter, dimension(6) :: components = [DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ]
    integer, parameter, dimension(2,2) :: connectivity = reshape([1, 2, 3, 4], [2,2])
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: release_angle = 20.0d0 * pi / 180.0d0
    real(real64), parameter :: length = 1.0d0
    real(real64), parameter :: radius = 0.05d0
    real(real64), parameter :: gravity = 9.81d0

    type :: constant_load
        real(real64), allocatable, dimension(:) :: force
    end type

    type(linear_elastic_isotropic_material), target :: material
    type(euler_bernoulli_beam_element) :: beam
    type(finite_element_mesh) :: mesh
    type(flexible_transient_settings) :: settings
    type(nodal_condition), dimension(12) :: supports
    type(degree_of_freedom), allocatable, dimension(:) :: root
    type(degree_of_freedom), allocatable, dimension(:) :: tip
    type(degree_of_freedom), allocatable, dimension(:) :: vertical
    type(multiplot) :: figure
    type(plot_2d) :: angle_plot
    type(plot_2d) :: bending_plot
    type(constant_load) :: gravity_load
    real(real64), dimension(3) :: axis
    real(real64), dimension(2) :: chord
    real(real64), dimension(nsamples) :: times
    real(real64), dimension(nsamples) :: angle
    real(real64), dimension(nsamples) :: bending
    real(real64), allocatable, dimension(:) :: initial
    real(real64), allocatable, dimension(:) :: gravity_acceleration
    real(real64), allocatable, dimension(:) :: constraints
    real(real64), allocatable, dimension(:,:) :: gradient
    real(real64), allocatable, dimension(:,:) :: mass
    real(real64), allocatable, dimension(:,:) :: q
    real(real64), allocatable, dimension(:,:) :: v
    real(real64) :: max_drift
    integer :: i
    integer :: ndof

    material%modulus = 1.0d6
    material%poissons_ratio = 0.3d0
    material%density = 850.0d0
    beam%material => material
    beam%area = pi * radius**2
    beam%iyy = pi * radius**4 / 4.0d0
    beam%izz = beam%iyy
    beam%j = 2.0d0 * beam%iyy

    ! Nodes 1 and 2 form a clamped support stub. Node 3 is the coincident
    ! pendulum root, joined to node 2 by a revolute joint about global Z.
    axis = [sin(release_angle), -cos(release_angle), 0.0d0]
    mesh = finite_element_mesh(4, 2)
    mesh%nodes = [fe_node(1, [-0.2d0, 0.0d0, 0.0d0]), &
        fe_node(2, [0.0d0, 0.0d0, 0.0d0]), fe_node(3, [0.0d0, 0.0d0, 0.0d0]), &
        fe_node(4, length * axis)]
    do i = 1, size(connectivity,2)
        beam%node_indices = connectivity(:,i)
        call mesh%assign_element(i, beam)
    end do
    mesh%joints = [flexible_joint(REVOLUTE_JOINT, 2, 3)]
    do i = 1, size(components)
        supports(i) = nodal_condition(1, components(i))
        supports(i+6) = nodal_condition(2, components(i))
    end do

    ! Resolve physical DOFs rather than assuming global equation numbers.
    call mesh%resolve_dofs([3, 3, 3], [DOF_UX, DOF_UY, DOF_RZ], root, ndof)
    call mesh%resolve_dofs([4, 4], [DOF_UX, DOF_UY], tip)
    call mesh%resolve_dofs([1, 2, 3, 4], [DOF_UY, DOF_UY, DOF_UY, DOF_UY], vertical)
    allocate(initial(ndof), gravity_acceleration(ndof), source=0.0d0)
    ! The release angle is already in the undeformed mesh. Convert uniform
    ! downward acceleration to a consistent nodal load with the reference mass.
    gravity_acceleration(vertical%equation_id) = -gravity
    call mesh%assemble_mass(mass)
    gravity_load%force = matmul(mass, gravity_acceleration)

    times = [(0.01d0 * (i-1), i=1,nsamples)]
    settings%solver = FLEXIBLE_KENNEDY_CARPENTER_5
    settings%absolute_tolerance = 1.0e-7
    settings%relative_tolerance = 1.0e-5
    call solve_flexible_transient(mesh, times, supports, initial, initial, q, v, settings, &
        load=apply_gravity, args=gravity_load)

    max_drift = 0.0d0
    angle = release_angle + q(root(3)%equation_id,:)
    do i = 1, nsamples
        chord = length * axis(:2) + q(tip%equation_id,i) - q(root(:2)%equation_id,i)
        ! Project the tip relative to the pivot onto the rotating root normal.
        bending(i) = dot_product(chord, [cos(angle(i)), sin(angle(i))])
        call mesh%assemble_constraints(q(:,i), constraints, gradient)
        max_drift = max(max_drift, maxval(abs(constraints)))
    end do
    print "(A,F8.3,A)", "Peak transverse tip deflection: ", 1000.0d0 * maxval(abs(bending)), " mm"
    print "(A,ES10.3)", "Maximum joint position residual: ", max_drift

    call figure%initialize(2, 1, width=1000, height=700)
    call angle_plot%initialize()
    call angle_plot%set_title("Flexible pendulum released from 20 degrees")
    call angle_plot%set_y_axis_title("Pivot angle (deg)")
    call angle_plot%push(times, angle * 180.0d0 / pi, lc=CLR_BLUE)
    call bending_plot%initialize()
    call bending_plot%set_title("Bending relative to the rotating root")
    call bending_plot%set_y_axis_title("Transverse tip deflection (mm)")
    call bending_plot%set_x_axis_title("Time (s)")
    call bending_plot%push(times, bending * 1000.0d0, lc=CLR_RED)
    call figure%set(1, 1, angle_plot)
    call figure%set(2, 1, bending_plot)
    call figure%draw()

contains
    subroutine apply_gravity(time, force, args)
        !! Return the precomputed full mesh load, independent of time.
        real(real64), intent(in) :: time
        real(real64), intent(out), dimension(:) :: force
        class(*), intent(in), optional :: args

        if (.not.present(args)) error stop "Missing pendulum gravity force"
        select type (args)
        type is (constant_load)
            force = args%force
        class default
            error stop "Invalid pendulum gravity force"
        end select
    end subroutine
end program
