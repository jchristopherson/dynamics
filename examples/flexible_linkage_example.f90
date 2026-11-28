program flexible_linkage_example
    !! Two gravity-driven Timoshenko links with a sprung and damped tip.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dynamics, only : fe_node, degree_of_freedom, finite_element_mesh, timoshenko_beam_element, &
        linear_elastic_isotropic_material, nodal_condition, flexible_joint, REVOLUTE_JOINT, &
        flexible_spring_element, flexible_damper_element, flexible_point_mass_element, &
        DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ, flexible_transient_settings, &
        solve_flexible_transient, FLEXIBLE_KENNEDY_CARPENTER_5
    use fplot_core, only : multiplot, plot_2d, legend, GNUPLOT_TERMINAL_PNG, CLR_BLUE, CLR_RED
    implicit none

    integer, parameter :: nsamples = 251
    integer, parameter, dimension(6) :: components = [DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ]
    integer, parameter, dimension(3) :: translations = [DOF_UX, DOF_UY, DOF_UZ]
    integer, parameter, dimension(2,3) :: beam_nodes = reshape([1, 2, 3, 4, 5, 6], [2,3])
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: base_release = 25.0d0*pi/180.0d0
    real(real64), parameter :: tip_release = -10.0d0*pi/180.0d0
    real(real64), parameter :: damping = 1.5d0

    type :: constant_load
        real(real64), allocatable, dimension(:) :: force
    end type

    type(linear_elastic_isotropic_material), target :: material
    type(timoshenko_beam_element) :: beam
    type(finite_element_mesh) :: mesh
    type(flexible_transient_settings) :: settings
    type(constant_load) :: gravity
    type(nodal_condition), dimension(15) :: supports
    type(degree_of_freedom), allocatable, dimension(:) :: hinges
    type(degree_of_freedom), allocatable, dimension(:) :: tip
    type(degree_of_freedom), allocatable, dimension(:) :: vertical
    type(multiplot) :: figure
    type(plot_2d) :: hinge_plot
    type(plot_2d) :: extension_plot
    class(legend), pointer :: hinge_legend
    real(real64), dimension(3,3) :: tip_inertia
    real(real64), dimension(3) :: elbow
    real(real64), dimension(3) :: chord
    real(real64), dimension(3) :: reference_chord
    real(real64), dimension(nsamples) :: times
    real(real64), dimension(nsamples) :: base_angle
    real(real64), dimension(nsamples) :: elbow_angle
    real(real64), dimension(nsamples) :: extension
    real(real64), allocatable, dimension(:) :: initial
    real(real64), allocatable, dimension(:) :: acceleration
    real(real64), allocatable, dimension(:) :: g
    real(real64), allocatable, dimension(:,:) :: jac
    real(real64), allocatable, dimension(:,:) :: mass
    real(real64), allocatable, dimension(:,:) :: q
    real(real64), allocatable, dimension(:,:) :: v
    real(real64) :: reference_length
    real(real64) :: extension_rate
    real(real64) :: peak_power
    real(real64) :: max_drift
    integer :: ndof
    integer :: i

    material = linear_elastic_isotropic_material(2.0d4, 0.3d0, 35.0d0)
    beam%material => material
    beam%area = 0.04d0
    beam%iyy = 2.0d-4
    beam%izz = 2.0d-4
    beam%j = 4.0d-4
    tip_inertia = 0.0d0
    do i = 1, 3
        tip_inertia(i,i) = 1.0d-3
    end do

    ! Coincident node pairs at the base and elbow are connected by joints,
    ! not merged: the two sides of each hinge can rotate independently.
    elbow = 0.7d0*[sin(base_release), -cos(base_release), 0.0d0]
    mesh = finite_element_mesh(7, 6)
    mesh%nodes = [fe_node(1, [-0.2d0, 0.0d0, 0.0d0]), fe_node(2, [0.0d0, 0.0d0, 0.0d0]), &
        fe_node(3, [0.0d0, 0.0d0, 0.0d0]), fe_node(4, elbow), fe_node(5, elbow), &
        fe_node(6, elbow + 0.55d0*[sin(tip_release), -cos(tip_release), 0.0d0]), &
        fe_node(7, [-0.65d0, -0.9d0, 0.0d0])]
    do i = 1, size(beam_nodes,2)
        beam%node_indices = beam_nodes(:,i)
        call mesh%assign_element(i, beam)
    end do
    call mesh%assign_element(4, flexible_spring_element(7, 6, 18.0d0))
    call mesh%assign_element(5, flexible_damper_element(7, 6, damping))
    call mesh%assign_element(6, flexible_point_mass_element(6, 0.25d0, tip_inertia))
    mesh%joints = [flexible_joint(REVOLUTE_JOINT, 2, 3), flexible_joint(REVOLUTE_JOINT, 4, 5)]
    do i = 1, size(components)
        supports(i) = nodal_condition(1, components(i))
        supports(i+6) = nodal_condition(2, components(i))
    end do
    do i = 1, size(translations)
        supports(i+12) = nodal_condition(7, translations(i))
    end do

    call mesh%resolve_dofs([3, 4, 5], [DOF_RZ, DOF_RZ, DOF_RZ], hinges, ndof)
    call mesh%resolve_dofs([6, 6, 6], translations, tip)
    call mesh%resolve_dofs([(i, i=1,7)], [(DOF_UY, i=1,7)], vertical)
    allocate(initial(ndof), acceleration(ndof), source=0.0d0)
    acceleration(vertical%equation_id) = -9.81d0
    call mesh%assemble_mass(mass)
    gravity%force = matmul(mass, acceleration)
    times = [(0.01d0*(i-1), i=1,nsamples)]
    settings%solver = FLEXIBLE_KENNEDY_CARPENTER_5
    settings%nthreads = 8
    call solve_flexible_transient(mesh, times, supports, initial, initial, q, v, settings, &
        load=apply_gravity, args=gravity)

    base_angle = (base_release + q(hinges(1)%equation_id,:))*180.0d0/pi
    elbow_angle = (tip_release - base_release + q(hinges(3)%equation_id,:) - &
        q(hinges(2)%equation_id,:))*180.0d0/pi
    reference_chord = mesh%nodes(6)%coordinates - mesh%nodes(7)%coordinates
    reference_length = norm2(reference_chord)
    max_drift = 0.0d0
    peak_power = 0.0d0
    do i = 1, nsamples
        chord = reference_chord + q(tip%equation_id,i)
        extension(i) = norm2(chord) - reference_length
        extension_rate = dot_product(chord, v(tip%equation_id,i))/norm2(chord)
        peak_power = max(peak_power, damping*extension_rate**2)
        call mesh%assemble_constraints(q(:,i), g, jac)
        max_drift = max(max_drift, maxval(abs(g)))
    end do
    if (.not.all(ieee_is_finite(extension))) error stop "Nonfinite linkage extension"
    print "(A,ES10.3)", "Maximum sampled joint position residual: ", max_drift
    print "(A,F8.3,A)", "Peak damper power dissipation: ", peak_power, " W"

    call figure%initialize(2, 1, width=1000, height=700)
    call hinge_plot%initialize()
    call hinge_plot%set_title("Gravity-driven Timoshenko linkage")
    call hinge_plot%set_y_axis_title("Hinge angle (deg)")
    hinge_legend => hinge_plot%get_legend()
    call hinge_legend%set_is_visible(.true.)
    call hinge_plot%push(times, base_angle, name="Base")
    call hinge_plot%push(times, elbow_angle, name="Elbow")
    call extension_plot%initialize()
    call extension_plot%set_title("Follower spring and damper")
    call extension_plot%set_y_axis_title("Connector extension (mm)")
    call extension_plot%set_x_axis_title("Time (s)")
    call extension_plot%push(times, 1000.0d0*extension)
    call figure%set(1, 1, hinge_plot)
    call figure%set(2, 1, extension_plot)
    call figure%draw()

contains
    subroutine apply_gravity(time, force, args)
        !! Return the preassembled constant gravity force.
        real(real64), intent(in) :: time
        real(real64), intent(out), dimension(:) :: force
        class(*), intent(in), optional :: args
        if (.not.present(args)) error stop "Missing linkage gravity load"
        select type (args)
        type is (constant_load)
            force = args%force
        class default
            error stop "Invalid linkage gravity load"
        end select
    end subroutine
end program
