! This example simulates a planar double pendulum whose rods are flexible.
! Each rod is an Euler-Bernoulli beam, discretized into a chain of spatial beam
! elements and assembled in the multibody framework:
!   - a spherical joint pins the proximal end of rod 1 to ground,
!   - a spherical joint connects the distal end of rod 1 to the proximal end of rod 2,
!   - gravity acts on the consistent element mass matrices.
! The rods start horizontal and at rest. The motion of the connecting joint and
! the free end is tracked and plotted with FPLOT.
program example
    use iso_fortran_env, only : int32, real64
    use dynamics
    use fplot_core
    implicit none

    ! Discretization and geometry
    integer(int32), parameter :: nel = 4
        ! Number of beam elements used to represent each rod.
    real(real64), parameter :: rod_length = 1.0d0
        ! Length of each rod [m].
    real(real64), parameter :: width = 5.0d-3
        ! Side length of the square rod cross section [m].

    ! Material properties (aluminum)
    real(real64), parameter :: youngs = 70.0d9
        ! Young's modulus [Pa].
    real(real64), parameter :: poisson = 0.33d0
        ! Poisson's ratio [-].
    real(real64), parameter :: density = 2700.0d0
        ! Mass density [kg/m^3].

    ! Environment and simulation controls
    real(real64), parameter :: gravity = 9.80665d0
        ! Gravitational acceleration magnitude [m/s^2].
    real(real64), parameter :: dt = 1.0d-3
        ! Fixed time step [s].
    integer(int32), parameter :: nstep = 10000
        ! Number of time steps; 10 seconds of motion.

    ! Derived quantities
    real(real64), parameter :: h = rod_length / nel
        ! Length of each beam element [m].
    real(real64), parameter :: area = width**2
        ! Cross-sectional area [m^2].
    real(real64), parameter :: inertia = width**4 / 12.0d0
        ! Second moment of area of a square section, the same about both axes [m^4].
    real(real64), parameter :: torsion = 0.1406d0 * width**4
        ! Saint-Venant torsion constant of a square section [m^4].
    integer(int32), parameter :: rod_dofs = 6 + 6*nel
        ! Coordinates per rod: six carrier coordinates plus six per beam node after the root.

    ! Local variables
    integer(int32) :: body
    integer(int32) :: e
    integer(int32) :: i
    integer(int32) :: ndof
    real(real64), allocatable, dimension(:) :: q
    real(real64), allocatable, dimension(:) :: v
    real(real64), allocatable, dimension(:) :: time
    real(real64), allocatable, dimension(:) :: x_joint, y_joint
    real(real64), allocatable, dimension(:) :: x_free, y_free
    real(real64), dimension(3) :: joint_position
    real(real64), dimension(3) :: free_position
    type(rigid_body) :: carrier
    type(flexible_euler_bernoulli_3d) :: beam
    type(dynamic_joint) :: link_joint
    type(multibody_model) :: model
    type(multibody_integrator) :: integrator
    type(multibody_state) :: state
    type(multibody_state), allocatable, dimension(:) :: history
    type(plot_2d) :: plt_path
    type(plot_2d) :: plt_time
    type(legend), pointer :: lgnd

    ! Build the model. Each rod is a carrier frame, which has zero rigid mass
    ! and inertia, with flexible beam elements attached. The carrier holds the
    ! rigid-body motion (origin plus rotation vector), and the elastic
    ! coordinates hold the bending. The mass lives entirely in the beam elements.
    carrier = rigid_body(m = 0.0d0)
    carrier%inertia = 0.0d0
    call model%add_body(carrier, 6*nel)
    call model%add_body(carrier, 6*nel)

    ! Section and material data shared by every beam element. The rods lie
    ! along the global x axis, so the default local y axis of the beam is valid.
    beam%material = linear_elastic(youngs, poisson, density)
    beam%area = area
    beam%second_moment_of_area = [inertia, inertia]
    beam%torsion_constant = torsion

    ! Attach a chain of nel elements to each rod. Element e spans nodes e-1
    ! to e. Adjacent elements share the node between them by reusing its
    ! elastic coordinates, which keeps each rod a single continuous beam.
    do body = 1, 2
        do e = 1, nel
            beam%reference_coordinates = 0.0d0
            beam%reference_coordinates(1,1) = (e-1) * h
            beam%reference_coordinates(1,2) = e * h
            call model%add_element(body, beam, [node_dofs(e-1), node_dofs(e)])
        end do
    end do

    ! Joint 1: spherical joint between ground (body 0) and the proximal end of
    ! rod 1. Both frames are identity transforms at the origin, and the rod
    ! root has no elastic motion, so no deformation maps are supplied.
    link_joint = make_joint(SPHERICAL_JOINT, 0, 1)
    call model%add_joint(link_joint)

    ! Joint 2: spherical joint at the connecting point. The rod 1 side is the
    ! distal end, at translation rod_length along its carrier x axis. Its
    ! deformation is read from the tip node, so the parent map is that node's
    ! elastic indices. The rod 2 side is its carrier origin, its proximal end.
    link_joint = make_joint(SPHERICAL_JOINT, 1, 2)
    link_joint%parent_frame(1,4) = rod_length
    call model%add_joint(link_joint, parent_map = node_dofs(nel))

    ! Gravity points along the negative global y axis.
    model%gravity = [0.0d0, -gravity, 0.0d0]

    ! Initial configuration: horizontal and at rest. Rod 1 lies along +x from
    ! the pivot at the origin, and rod 2 starts at the tip of rod 1, also along
    ! +x. The rotation vectors are zero and the beams are undeformed.
    ndof = model%get_dof_count()
    allocate(q(ndof), v(ndof))
    q = 0.0d0
    v = 0.0d0
    q(rod_dofs+1:rod_dofs+3) = [rod_length, 0.0d0, 0.0d0]

    ! Initialize the state, which checks the constraints and computes consistent
    ! accelerations, then integrate for nstep steps. The returned history
    ! holds nstep+1 samples, starting with the initial state.
    call integrator%initialize(model, state, 0.0d0, q, v)
    history = integrator%solve(model, state, dt, nstep)

    ! Recover the global positions of the connecting joint (distal end of rod 1)
    ! and of the free end (distal end of rod 2) at each sample.
    allocate(time(nstep+1), x_joint(nstep+1), y_joint(nstep+1))
    allocate(x_free(nstep+1), y_free(nstep+1))
    do i = 1, nstep+1
        time(i) = history(i)%time
        joint_position = rod_tip(history(i)%position, 1)
        free_position = rod_tip(history(i)%position, 2)
        x_joint(i) = joint_position(1)
        y_joint(i) = joint_position(2)
        x_free(i) = free_position(1)
        y_free(i) = free_position(2)
    end do

    print "(A,2F10.4)", "Final connecting joint position [m]: ", &
        x_joint(nstep+1), y_joint(nstep+1)
    print "(A,2F10.4)", "Final free end position [m]: ", &
        x_free(nstep+1), y_free(nstep+1)

    ! Plot the planar path of both points.
    call plt_path%initialize()
    call plt_path%set_x_axis_title("x [m]")
    call plt_path%set_y_axis_title("y [m]")
    lgnd => plt_path%get_legend()
    call lgnd%set_is_visible(.true.)
    call plt_path%push(x_joint, y_joint, name = "Connecting joint")
    call plt_path%push(x_free, y_free, name = "Free end")
    call plt_path%draw()

    ! Plot the vertical position of both points over time.
    call plt_time%initialize()
    call plt_time%set_x_axis_title("Time [s]")
    call plt_time%set_y_axis_title("Vertical Position [m]")
    lgnd => plt_time%get_legend()
    call lgnd%set_is_visible(.true.)
    call plt_time%push(time, y_joint, name = "Connecting joint")
    call plt_time%push(time, y_free, name = "Free end")
    call plt_time%draw()

contains
! ------------------------------------------------------------------------------
    pure function node_dofs(node) result(map)
        ! Maps the six coordinates of a beam node to the rod's elastic indices.
        ! Node zero is the rod root, welded to its carrier, so a zero entry
        ! restrains all six of its coordinates.
        integer(int32), intent(in) :: node
            ! Beam node index, from zero (root) to nel (distal end).
        integer, dimension(6) :: map
            ! Elastic indices for the three translations and three rotations.

        if (node == 0) then
            map = 0
        else
            map = 6*(node-1) + [1, 2, 3, 4, 5, 6]
        end if
    end function

! ------------------------------------------------------------------------------
    function make_joint(kind, parent, child) result(joint)
        ! Creates a joint whose attachment frames are identity transforms.
        ! Frame translations are in column four of each 4-by-4 transform.
        integer(int32), intent(in) :: kind
            ! Joint type constant, such as SPHERICAL_JOINT.
        integer(int32), intent(in) :: parent
            ! Parent body index, or zero for ground.
        integer(int32), intent(in) :: child
            ! Child body index.
        type(dynamic_joint) :: joint
            ! The joint descriptor.
        integer :: i

        joint%joint_type = kind
        joint%parent_body = parent
        joint%child_body = child
        joint%parent_frame = 0.0d0
        joint%child_frame = 0.0d0
        do i = 1, 4
            joint%parent_frame(i,i) = 1.0d0
            joint%child_frame(i,i) = 1.0d0
        end do
    end function

! ------------------------------------------------------------------------------
    function rod_tip(q, body) result(p)
        ! Returns the global position of the distal end of a rod. This is the
        ! carrier transform applied to the tip node's reference position plus
        ! its elastic translation, the same kinematics the multibody model uses.
        real(real64), intent(in), dimension(:) :: q
            ! Assembled coordinates at one sample.
        integer(int32), intent(in) :: body
            ! Rod index: 1 for the first rod, 2 for the second.
        real(real64), dimension(3) :: p
            ! Global position of the distal end.
        integer(int32) :: lo
            ! Index of the first coordinate of this rod's block.
        real(real64), dimension(3) :: u
            ! Elastic translation of the tip node, in the carrier frame.
        real(real64), dimension(3,3) :: r
            ! Carrier rotation matrix.

        lo = (body-1)*rod_dofs + 1
        ! The tip node's elastic coordinates are the last six of the block.
        ! Their first three entries are translations.
        u = q(lo+6*nel+1:lo+6*nel+3)
        r = carrier_rotation(q(lo+3:lo+5))
        p = q(lo:lo+2) + matmul(r, [rod_length, 0.0d0, 0.0d0] + u)
    end function

! ------------------------------------------------------------------------------
    pure function carrier_rotation(p) result(r)
        ! Rotation matrix for a rotation vector p, computed with the exponential
        ! map (Rodrigues' formula). This matches the convention the multibody
        ! solver uses for carrier rotations.
        real(real64), intent(in), dimension(3) :: p
            ! Rotation vector in radians.
        real(real64), dimension(3,3) :: r
            ! The corresponding rotation matrix.
        real(real64), dimension(3,3) :: k
            ! Cross-product matrix of p.
        real(real64), dimension(3,3) :: eye
            ! Identity matrix.
        real(real64) :: angle
            ! Rotation angle, the magnitude of p.
        real(real64) :: s
            ! Coefficient of the first-order term, sin(angle)/angle.
        real(real64) :: c
            ! Coefficient of the second-order term, (1-cos(angle))/angle^2.

        eye = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3,3])
        k = reshape([0.0d0, p(3), -p(2), -p(3), 0.0d0, p(1), p(2), -p(1), 0.0d0], [3,3])
        angle = norm2(p)
        if (angle < 1.0d-8) then
            ! Series limits avoid dividing by a vanishing angle.
            s = 1.0d0
            c = 0.5d0
        else
            s = sin(angle) / angle
            c = (1.0d0 - cos(angle)) / angle**2
        end if
        r = eye + s*k + c*matmul(k,k)
    end function

end program