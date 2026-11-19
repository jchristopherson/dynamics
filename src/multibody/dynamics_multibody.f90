! SPDX-License-Identifier: MIT
module dynamics_multibody
    !! Spatial floating-frame assembly for rigid bodies and flexible elements.
    !! Each body's coordinates are `[origin, rotation_vector, elastic_dofs]`.
    !! Origins are in world coordinates; elastic coordinates are body-local.
    !! The carrier rotation vector is expressed in world axes. Its time
    !! derivative is a coordinate rate, not angular velocity.
    !! Element inertia uses constant reference mass pulled back through the
    !! current nodal velocity maps; for deformed beams this is an approximation.
    use iso_fortran_env, only : real64, int32
    use ieee_arithmetic, only : ieee_is_finite
    use dynamics_rigid_bodies, only : rigid_body
    use dynamics_flexible_types, only : flexible_element, finite_element, dissipative_element
    use dynamics_flexible_truss, only : flexible_truss_3d
    use dynamics_flexible_discrete, only : flexible_spring_3d
    use dynamics_flexible_beams, only : flexible_euler_bernoulli_3d, flexible_timoshenko_3d, &
        flexible_euler_bernoulli_2d, flexible_timoshenko_2d
    use dynamics_linkage_dynamics, only : dynamic_joint, joint_constraint_count, append_spatial_joint_constraints
    use dynamics_kinematics, only : transform_inverse
    implicit none
    private
    public :: multibody_model, multibody_constraint, multibody_load

    abstract interface
        subroutine multibody_constraint(time, q, value)
            !! Smooth additional holonomic residuals, equal to zero on the constraint manifold.
            import real64
            real(real64), intent(in) :: time
                !! Evaluation time, including for prescribed motion.
            real(real64), intent(in), dimension(:) :: q
                !! Current coordinates in model order.
            real(real64), intent(out), dimension(:) :: value
                !! All custom residual rows, in the order used for their multipliers.
        end subroutine

        subroutine multibody_load(time, q, v, force)
            !! Applied generalized forces conjugate to the assembled coordinate rates.
            import real64
            real(real64), intent(in) :: time
                !! Load evaluation time.
            real(real64), intent(in), dimension(:) :: q
                !! Current coordinates in model order.
            real(real64), intent(in), dimension(:) :: v
                !! Current coordinate rates in model order.
            real(real64), intent(out), dimension(:) :: force
                !! One generalized force per coordinate, excluding model gravity.
        end subroutine
    end interface

    type :: element_entry
        !! An owned element and its cached reference data for assembly.
        class(flexible_element), allocatable :: element
            !! Element copy preserving its concrete type and constitutive data.
        integer, allocatable, dimension(:) :: map
            !! Element DOF to body-local elastic DOF map; zero denotes restraint.
        real(real64), allocatable, dimension(:,:) :: mass
            !! Reference mass matrix in element DOF order.
        real(real64), dimension(3,2) :: reference = 0.0d0
            !! Body-local reference node positions, one per column.
    end type

    type :: body_entry
        !! One carrier and its optional elastic coordinate block.
        type(rigid_body) :: carrier
            !! Mass properties of the carrier alone, excluding element mass.
        integer :: first
            !! Global index of the first carrier translation.
        integer :: elastic_count
            !! Number of active elastic coordinates following the six carrier coordinates.
        type(element_entry), allocatable, dimension(:) :: elements
            !! Elements expressed in this carrier's frame.
    end type

    type :: joint_entry
        !! A spatial joint with optional deformation at either attachment.
        type(dynamic_joint) :: joint
            !! Joint kind, endpoint bodies and undeformed attachment frames.
        integer, dimension(6) :: parent_map
            !! Parent elastic indices in translation/rotation order.
        integer, dimension(6) :: child_map
            !! Child elastic indices in translation/rotation order.
    end type

    type :: multibody_model
        !! Assembles a spatial floating-frame system without imposing a time integrator.
        !! Bodies and joint rows are stored in insertion order. Ground has index zero.
        !! The coordinate count includes constrained coordinates; joint and custom
        !! residuals enforce restrictions without eliminating these coordinates.
        private
        type(body_entry), allocatable, dimension(:) :: bodies
            !! Moving bodies in coordinate-block order.
        type(joint_entry), allocatable, dimension(:) :: joints
            !! Joint descriptors in residual-row order.
        integer :: ndof = 0
            !! Number of assembled coordinates.
        integer :: njoint_rows = 0
            !! Number of scalar rows supplied by spatial joints.
        integer :: ncustom = 0
            !! Number of custom rows following the joint rows.
        procedure(multibody_constraint), pointer, nopass :: custom => null()
            !! Optional callback filling the custom constraint rows.
        procedure(multibody_load), pointer, nopass :: load => null()
            !! Optional callback filling the applied generalized load vector.
        real(real64), public, dimension(3) :: gravity = 0.0d0
            !! Uniform world-frame acceleration, applied to carrier and element masses.
    contains
        procedure :: add_body
        procedure :: add_element
        procedure :: add_joint
        procedure :: set_constraints
        procedure :: set_load
        procedure :: get_dof_count
        procedure :: get_constraint_count
        procedure :: get_body_coordinates
        procedure :: mass_matrix
        procedure :: equations
        procedure :: constraints
    end type

contains

subroutine add_body(this, carrier, elastic_dofs)
    !! Appends a carrier with optional elastic DOFs. For a flexible-only body,
    !! pass zero carrier mass and inertia to avoid double counting element mass.
    class(multibody_model), intent(inout) :: this
        !! Model receiving the new body at the end of its coordinate layout.
    type(rigid_body), intent(in) :: carrier
        !! Carrier mass, center of mass in local coordinates, and inertia about that center.
    integer, intent(in), optional :: elastic_dofs
        !! Number of active body-local elastic coordinates; defaults to zero.
    type(body_entry) :: entry

    entry%carrier = carrier
    entry%first = this%ndof + 1
    entry%elastic_count = 0
    if (present(elastic_dofs)) entry%elastic_count = elastic_dofs
    if (entry%elastic_count < 0 .or. .not. ieee_is_finite(carrier%mass) .or. carrier%mass < 0.0d0 .or. &
        .not. all(ieee_is_finite(carrier%cg)) .or. .not. all(ieee_is_finite(carrier%inertia))) &
        error stop "Invalid multibody carrier or elastic DOF count"
    if (maxval(abs(carrier%inertia-transpose(carrier%inertia))) > &
        100*epsilon(1.0d0)*maxval(abs(carrier%inertia))) error stop "Carrier inertia must be symmetric"
    if (.not. allocated(this%bodies)) allocate(this%bodies(0))
    this%bodies = [this%bodies, entry]
    this%ndof = this%ndof + 6 + entry%elastic_count
end subroutine

subroutine add_element(this, body, element, map, reference)
    !! Copies a spatial two-node element into a floating frame. `map` contains
    !! body-local elastic indices in element order; zero restrains that DOF.
    !! `reference` supplies the two reference positions (columns) of custom
    !! massive elements. The built-in spatial elements provide this geometry.
    class(multibody_model), intent(inout) :: this
        !! Model that owns the specified body.
    integer, intent(in) :: body
        !! One-based index of the moving body receiving the element.
    integer, intent(in), dimension(:) :: map
        !! One entry per element DOF. Reuse an index to share an active coordinate.
    class(flexible_element), intent(in) :: element
        !! Spatial two-node element to copy, including its material and geometry.
    real(real64), intent(in), optional, dimension(3,2) :: reference
        !! Reference positions for a custom massive element, in the carrier frame.
    type(element_entry) :: entry
    integer :: n

    call check_body(this, body)
    n = element%get_dof_count()
    if ((n /= 6 .and. n /= 12) .or. size(map) /= n) error stop "Expected a spatial two-node element"
    call check_map(map, this%bodies(body)%elastic_count)
    allocate(entry%element, source=element)
    entry%map = map
    entry%mass = element%mass_matrix()
    if (any(shape(entry%mass) /= [n,n]) .or. .not. all(ieee_is_finite(entry%mass))) &
        error stop "Invalid flexible element mass matrix"
    if (maxval(abs(entry%mass-transpose(entry%mass))) > &
        100*epsilon(1.0d0)*maxval(abs(entry%mass))) error stop "Element mass matrix must be symmetric"
    select type (element)
    class is (flexible_truss_3d)
        entry%reference = element%reference_coordinates
    class is (flexible_spring_3d)
        entry%reference = element%reference_coordinates
    class is (flexible_euler_bernoulli_3d)
        entry%reference = element%reference_coordinates
    class is (flexible_timoshenko_3d)
        entry%reference = element%reference_coordinates
    class is (flexible_euler_bernoulli_2d)
        error stop "Planar beams cannot be added to a spatial multibody model"
    class is (flexible_timoshenko_2d)
        error stop "Planar beams cannot be added to a spatial multibody model"
    class default
        if (any(abs(entry%mass) > 0.0d0) .and. .not. present(reference)) &
            error stop "Custom massive elements require reference coordinates"
        if (present(reference)) entry%reference = reference
    end select
    if (.not. all(ieee_is_finite(entry%reference))) error stop "Invalid element reference coordinates"
    associate(b => this%bodies(body))
        if (.not. allocated(b%elements)) allocate(b%elements(0))
        b%elements = [b%elements, entry]
    end associate
end subroutine

subroutine add_joint(this, joint, parent_map, child_map)
    !! Adds an existing spatial joint descriptor. For an elastic attachment,
    !! supply six body-local elastic indices `[ux,uy,uz,rx,ry,rz]` per end.
    !! A zero index fixes that component relative to the carrier.
    !! Each joint-frame origin is the attachment point before deformation.
    class(multibody_model), intent(inout) :: this
        !! Model receiving the joint after its existing joints.
    type(dynamic_joint), intent(in) :: joint
        !! Joint descriptor with proper homogeneous transforms and valid body indices.
    integer, intent(in), optional, dimension(6) :: parent_map
        !! Parent deformation map; omit for a rigid attachment or ground.
    integer, intent(in), optional, dimension(6) :: child_map
        !! Child deformation map; omit for a rigid attachment or ground.
    type(joint_entry) :: entry

    if (joint%parent_body == joint%child_body) error stop "Joint must connect distinct bodies"
    entry%joint = joint
    entry%parent_map = 0
    entry%child_map = 0
    if (present(parent_map)) entry%parent_map = parent_map
    if (present(child_map)) entry%child_map = child_map
    call check_endpoint(joint%parent_body, entry%parent_map)
    call check_endpoint(joint%child_body, entry%child_map)
    call check_frame(joint%parent_frame)
    call check_frame(joint%child_frame)
    if (.not. allocated(this%joints)) allocate(this%joints(0))
    this%njoint_rows = this%njoint_rows + joint_constraint_count(joint%joint_type, .false.)
    this%joints = [this%joints, entry]
contains
    subroutine check_endpoint(body, map)
        !! Ground has no elastic motion; moving endpoints use their own coordinate map.
        integer, intent(in) :: body
            !! Endpoint body index, or zero for ground.
        integer, intent(in), dimension(6) :: map
            !! Elastic coordinate map for that endpoint.
        if (body == 0) then
            call check_map(map, 0)
        else
            call check_body(this, body)
            call check_map(map, this%bodies(body)%elastic_count)
        end if
    end subroutine
end subroutine

subroutine check_frame(frame)
    !! A spatial joint frame must be a proper rigid transform.
    use dynamics_helper, only : cross_product
    real(real64), intent(in), dimension(4,4) :: frame
        !! Transform to validate: proper rotation, finite translation and homogeneous last row.
    real(real64), dimension(3,3) :: r
    real(real64), dimension(3,3) :: eye
    real(real64), parameter :: tol = 1.0d-10

    r = frame(:3,:3)
    eye = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    if (.not. all(ieee_is_finite(frame))) error stop "Nonfinite multibody joint frame"
    if (maxval(abs(matmul(transpose(r),r)-eye)) > tol .or. &
        maxval(abs(cross_product(r(:,1),r(:,2))-r(:,3))) > tol .or. &
        maxval(abs(frame(4,:)-[0.0d0,0.0d0,0.0d0,1.0d0])) > tol) error stop "Invalid multibody joint frame"
end subroutine

subroutine set_constraints(this, count, evaluate)
    !! Registers additional, possibly time-dependent, constraint rows after
    !! all joint rows. With count zero and no callback, removes these rows.
    class(multibody_model), intent(inout) :: this
        !! Model whose custom constraints are replaced or cleared.
    integer, intent(in) :: count
        !! Number of residual rows to allocate for the callback.
    procedure(multibody_constraint), optional :: evaluate
        !! Callback required for positive `count`; its target must outlive model use.

    if (count < 0 .or. ((count > 0) .neqv. present(evaluate))) error stop "Invalid custom constraints"
    this%ncustom = count
    nullify(this%custom)
    if (present(evaluate)) this%custom => evaluate
end subroutine

subroutine set_load(this, evaluate)
    !! Registers a generalized load callback, or clears it if omitted.
    class(multibody_model), intent(inout) :: this
        !! Model whose applied-load callback is replaced or cleared.
    procedure(multibody_load), optional :: evaluate
        !! Generalized load callback; its target must outlive model use.

    nullify(this%load)
    if (present(evaluate)) this%load => evaluate
end subroutine

pure function get_dof_count(this) result(n)
    !! Number of assembled coordinates, before applying joint constraints.
    class(multibody_model), intent(in) :: this
        !! Model whose coordinate layout is queried.
    integer :: n
        !! Number of carrier and active elastic coordinates.
    n = this%ndof
end function

pure function get_constraint_count(this) result(n)
    !! Total number of scalar joint and user constraint rows.
    class(multibody_model), intent(in) :: this
        !! Model whose residual layout is queried.
    integer :: n
        !! Number of scalar constraint rows, including custom rows.
    n = this%njoint_rows + this%ncustom
end function

function get_body_coordinates(this, body) result(indices)
    !! Returns a body's global coordinate indices, including its carrier.
    class(multibody_model), intent(in) :: this
        !! Model providing the global coordinate layout.
    integer, intent(in) :: body
        !! One-based moving body index.
    integer, allocatable, dimension(:) :: indices
        !! Contiguous global indices: translation, rotation, then elastic coordinates.
    integer :: i
    integer :: first

    call check_body(this, body)
    first = this%bodies(body)%first
    indices = [(first+i-1, i=1,6+this%bodies(body)%elastic_count)]
end function

function constraints(this, time, q) result(value)
    !! Evaluates joint and user constraint residuals at the specified pose.
    class(multibody_model), intent(in) :: this
        !! Model providing the joints and optional callback.
    real(real64), intent(in) :: time
        !! Time passed to the custom constraint callback.
    real(real64), intent(in), dimension(:) :: q
        !! Assembled coordinates of length `get_dof_count()`.
    real(real64), dimension(this%get_constraint_count()) :: value
        !! Residuals in joint insertion order, followed by custom rows.
    real(real64), dimension(4,4) :: parent
    real(real64), dimension(4,4) :: child
    integer(int32) :: row
    integer :: i

    if (size(q) /= this%ndof) error stop "Wrong multibody coordinate count"
    row = 1
    if (allocated(this%joints)) then
        do i = 1, size(this%joints)
            associate(entry => this%joints(i))
                parent = attachment(this, q, entry%joint%parent_body, entry%joint%parent_frame, entry%parent_map)
                child = attachment(this, q, entry%joint%child_body, entry%joint%child_frame, entry%child_map)
                call append_spatial_joint_constraints(entry%joint%joint_type, &
                    matmul(transform_inverse(parent), child), value, row, use_half_angle=.true.)
            end associate
        end do
    end if
    if (associated(this%custom)) call this%custom(time, q, value(row:))
end function

function attachment(this, q, body, frame, map) result(world)
    !! Applies an elastic displacement/rotation before the carrier transform.
    class(multibody_model), intent(in) :: this
        !! Model providing the endpoint coordinate block.
    real(real64), intent(in), dimension(:) :: q
        !! Assembled coordinates.
    real(real64), intent(in), dimension(4,4) :: frame
        !! Undeformed carrier-local attachment, or world transform for ground.
    integer, intent(in) :: body
        !! Endpoint body index, or zero for ground.
    integer, intent(in), dimension(6) :: map
        !! Body-local elastic indices of the three translations and three rotations.
    real(real64), dimension(4,4) :: world
        !! Current world transform of the attachment.
    real(real64), dimension(3,3) :: r
    real(real64), dimension(3,3) :: jac
    real(real64), dimension(3,3) :: rd
    real(real64), dimension(6) :: u
    integer :: lo
    integer :: hi

    world = frame
    if (body == 0) return
    lo = this%bodies(body)%first
    hi = lo + 5 + this%bodies(body)%elastic_count
    u = gather(q(lo+6:hi), map)
    call rotation_kinematics(q(lo+3:lo+5), r, jac)
    call rotation_kinematics(u(4:6), rd, jac)
    world(1:3,1:3) = matmul(matmul(r, rd), frame(1:3,1:3))
    world(1:3,4) = q(lo:lo+2) + matmul(r, frame(1:3,4) + u(1:3))
end function

function mass_matrix(this, q) result(m)
    !! Returns the dense assembled metric for kinetic energy `0.5*v^T*M*v`.
    class(multibody_model), intent(in) :: this
        !! Model whose inertia is assembled.
    real(real64), intent(in), dimension(:) :: q
        !! Assembled coordinates at which to evaluate the velocity-to-motion maps.
    real(real64), dimension(this%ndof,this%ndof) :: m
        !! Dense mass matrix in assembled coordinate order, before applying joint constraints.
    integer :: i
    integer :: lo
    integer :: hi

    if (size(q) /= this%ndof) error stop "Wrong multibody coordinate count"
    m = 0.0d0
    if (.not. allocated(this%bodies)) return
    do i = 1, size(this%bodies)
        lo = this%bodies(i)%first
        hi = lo + 5 + this%bodies(i)%elastic_count
        m(lo:hi,lo:hi) = body_mass(this%bodies(i), q(lo:hi))
    end do
end function

function body_mass(body, q) result(m)
    !! Pulls reference nodal inertia back through instantaneous nodal velocity maps.
    type(body_entry), intent(in) :: body
        !! Carrier and its attached reference element data.
    real(real64), intent(in), dimension(:) :: q
        !! Coordinates of this body only, including its carrier.
    real(real64), dimension(size(q),size(q)) :: m
        !! Mass block including carrier-elastic coupling.
    real(real64), dimension(12,size(q)) :: b
    real(real64), dimension(3,3) :: r
    real(real64), dimension(3,3) :: jbody
    real(real64), dimension(3,3) :: jnode
    real(real64), dimension(3,3) :: rn
    real(real64), dimension(12) :: u
    real(real64), dimension(3) :: p
    real(real64), dimension(3,size(q)) :: linear
    integer :: i
    integer :: node
    integer :: n
    integer :: stride
    integer :: k
    integer :: a
    integer :: g

    call rotation_kinematics(q(4:6), r, jbody)
    jbody = matmul(transpose(r), jbody)
    linear = 0.0d0
    linear(:,1:3) = transpose(r)
    linear(:,4:6) = -matmul(skew(body%carrier%cg), jbody)
    m = body%carrier%mass * matmul(transpose(linear), linear)
    m(4:6,4:6) = m(4:6,4:6) + matmul(transpose(jbody), matmul(body%carrier%inertia, jbody))
    if (.not. allocated(body%elements)) return
    do i = 1, size(body%elements)
        associate(e => body%elements(i))
            if (.not. any(abs(e%mass) > 0.0d0)) cycle
            n = size(e%map)
            stride = n / 2
            u(:n) = gather(q(7:), e%map)
            b(:n,:) = 0.0d0
            do node = 1, 2
                k = (node-1)*stride
                p = e%reference(:,node) + u(k+1:k+3)
                b(k+1:k+3,1:3) = transpose(r)
                b(k+1:k+3,4:6) = -matmul(skew(p), jbody)
                do a = 1, 3
                    g = e%map(k+a)
                    if (g > 0) b(k+a,6+g) = 1.0d0
                end do
                if (stride == 6) then
                    call rotation_kinematics(u(k+4:k+6), rn, jnode)
                    b(k+4:k+6,4:6) = jbody
                    do a = 1, 3
                        g = e%map(k+3+a)
                        if (g > 0) b(k+4:k+6,6+g) = b(k+4:k+6,6+g) + jnode(:,a)
                    end do
                end if
            end do
            m = m + matmul(transpose(b(:n,:)), matmul(e%mass, b(:n,:)))
        end associate
    end do
end function

subroutine equations(this, time, q, v, m, force)
    !! Returns M and the right-hand side of `M*a + G^T*lambda = force`.
    !! The right-hand side includes gravity, loads, elastic/damping forces,
    !! and convective inertia `Mdot*v - grad_q(0.5*v^T*M*v)`.
    class(multibody_model), intent(in) :: this
        !! Model providing mass properties, gravity, elements and applied loads.
    real(real64), intent(in) :: time
        !! Time passed to the applied-load callback.
    real(real64), intent(in), dimension(:) :: q
        !! Coordinates in assembled order, of length `get_dof_count()`.
    real(real64), intent(in), dimension(:) :: v
        !! Coordinate rates in the same order as `q`.
    real(real64), intent(out), dimension(:,:) :: m
        !! Square mass matrix with one row and column per assembled coordinate.
    real(real64), intent(out), dimension(:) :: force
        !! Net generalized right-hand side, excluding joint reactions.
    integer :: i
    integer :: lo
    integer :: hi

    if (size(q) /= this%ndof .or. size(v) /= this%ndof .or. size(force) /= this%ndof .or. &
        any(shape(m) /= this%ndof)) error stop "Wrong multibody equation dimensions"
    m = 0.0d0
    force = 0.0d0
    if (associated(this%load)) call this%load(time, q, v, force)
    if (.not. allocated(this%bodies)) return
    do i = 1, size(this%bodies)
        lo = this%bodies(i)%first
        hi = lo + 5 + this%bodies(i)%elastic_count
        m(lo:hi,lo:hi) = body_mass(this%bodies(i), q(lo:hi))
        force(lo:hi) = force(lo:hi) + matmul(m(lo:hi,lo:lo+2), this%gravity) - &
            body_resistance(this%bodies(i), q(lo:hi), v(lo:hi))
    end do
end subroutine

function body_resistance(body, q, v) result(force)
    !! Evaluates the kinetic-energy derivatives by centered differences, then
    !! scatters the autodifferentiated elastic and viscous forces.
    type(body_entry), intent(in) :: body
        !! Carrier and its attached elastic and dissipative elements.
    real(real64), intent(in), dimension(:) :: q
        !! Carrier and elastic coordinates for this body only.
    real(real64), intent(in), dimension(:) :: v
        !! Coordinate rates for this body in the same order as `q`.
    real(real64), dimension(size(q)) :: force
        !! Convective inertia plus elastic and damping resistance, subtracted from applied loads.
    real(real64), dimension(size(q)) :: qp
    real(real64), dimension(size(q)) :: qm
    real(real64), dimension(size(v)) :: direction
    real(real64) :: h
    real(real64) :: speed
    real(real64), dimension(12) :: u
    real(real64), dimension(12) :: f
    real(real64), dimension(12,12) :: tangent
    integer :: i
    integer :: k
    integer :: n
    integer :: g

    force = 0.0d0
    direction = v
    direction(1:3) = 0.0d0   ! The mass metric is invariant under carrier translations.
    speed = maxval(abs(direction))
    if (speed > tiny(speed)) then
        h = epsilon(h)**(1.0d0/3.0d0)
        direction = direction / speed
        force = (speed/(2*h)) * matmul(body_mass(body, q+h*direction)-body_mass(body, q-h*direction), v)
    end if
    if (maxval(abs(v)) > 0.0d0) then
        do i = 4, size(q)
            h = epsilon(h)**(1.0d0/3.0d0)*max(1.0d0,abs(q(i)))
            qp = q
            qm = q
            qp(i) = qp(i) + h
            qm(i) = qm(i) - h
            force(i) = force(i) - dot_product(v, matmul(body_mass(body,qp)-body_mass(body,qm),v))/(4*h)
        end do
    end if
    if (.not. allocated(body%elements)) return
    do i = 1, size(body%elements)
        associate(e => body%elements(i))
            n = size(e%map)
            u(:n) = gather(q(7:), e%map)
            select type (element => e%element)
            class is (finite_element)
                call element%internal_force(u(:n), f(:n), tangent(:n,:n))
            class is (dissipative_element)
                call element%damping_matrix(u(:n), tangent(:n,:n))
                f(:n) = matmul(tangent(:n,:n), gather(v(7:), e%map))
            class default
                cycle
            end select
            do k = 1, n
                g = e%map(k)
                if (g > 0) force(6+g) = force(6+g) + f(k)
            end do
        end associate
    end do
end function

pure function gather(x, map) result(values)
    !! Gathers active elastic coordinates and inserts zeros at restrained DOFs.
    real(real64), intent(in), dimension(:) :: x
        !! Active body-local elastic coordinates or rates.
    integer, intent(in), dimension(:) :: map
        !! One-based indices into `x`, or zero for restrained components.
    real(real64), dimension(size(map)) :: values
        !! Values in map order, with zeros for restrained components.
    integer :: i

    values = 0.0d0
    do i = 1, size(map)
        if (map(i) > 0) values(i) = x(map(i))
    end do
end function

pure subroutine rotation_kinematics(p, r, left_jacobian)
    !! SO(3) exponential and left Jacobian: omega = J(p)*p_dot. The
    !! series avoids cancellation for small rotations.
    real(real64), intent(in), dimension(3) :: p
        !! Rotation vector in radians.
    real(real64), intent(out), dimension(3,3) :: r
        !! Rotation matrix obtained from the exponential map.
    real(real64), intent(out), dimension(3,3) :: left_jacobian
        !! Maps vector rates to angular velocity in the rotation's destination frame.
    real(real64) :: a
    real(real64) :: a2
    real(real64) :: s
    real(real64) :: c
    real(real64) :: d
    real(real64), dimension(3,3) :: k
    real(real64), dimension(3,3) :: k2
    real(real64), dimension(3,3) :: eye

    a2 = dot_product(p,p)
    if (a2 < 1.0d-6) then
        s = 1 - a2/6 + a2*a2/120
        c = 0.5d0 - a2/24 + a2*a2/720
        d = 1.0d0/6 - a2/120 + a2*a2/5040
    else
        a = sqrt(a2)
        s = sin(a)/a
        c = 2*(sin(a/2)/a)**2
        d = (1-s)/a2
    end if
    k = skew(p)
    k2 = matmul(k,k)
    eye = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    r = eye + s*k + c*k2
    left_jacobian = eye + c*k + d*k2
end subroutine

pure function skew(x) result(k)
    !! Cross-product matrix, such that `matmul(k,y)` is `x` cross `y`.
    real(real64), intent(in), dimension(3) :: x
        !! Left-hand vector of the cross product.
    real(real64), dimension(3,3) :: k
        !! Skew-symmetric matrix representing left multiplication by `x`.
    k = reshape([0.0d0,x(3),-x(2),-x(3),0.0d0,x(1),x(2),-x(1),0.0d0],[3,3])
end function

subroutine check_map(map, n)
    !! Validates a map before it is retained by the model.
    integer, intent(in), dimension(:) :: map
        !! Elastic indices to validate, with zero permitted.
    integer, intent(in) :: n
        !! Number of active coordinates available to the map.
    if (any(map < 0) .or. any(map > n)) error stop "Invalid multibody elastic coordinate map"
end subroutine

subroutine check_body(this, body)
    !! Validates a moving body index (ground is handled by the caller).
    class(multibody_model), intent(in) :: this
        !! Model whose moving body list is checked.
    integer, intent(in) :: body
        !! Expected one-based moving body index.
    if (.not. allocated(this%bodies)) error stop "No multibody carriers have been added"
    if (body < 1 .or. body > size(this%bodies)) error stop "Invalid moving body index"
end subroutine

end module
