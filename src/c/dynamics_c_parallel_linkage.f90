module dynamics_c_parallel_linkage
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use diffeq
    use spectrum, only : window
    use dynamics_error_handling
    use nonlin
    use dynamics_c_types
    implicit none

contains

function get_mechanism(obj) result(rst)
    ! Resolves an opaque mechanism handle.
    type(c_ptr), intent(in), value :: obj
    class(kinematic_mechanism), pointer :: rst

    type(c_mechanism_container), pointer :: cont

    rst => null()
    if (.not.c_associated(obj)) return
    call c_f_pointer(obj, cont)
    if (allocated(cont%item)) rst => cont%item
end function

! ------------------------------------------------------------------------------
function build_mechanism(planar, nlinks, links, njoints, joints, base, &
    effector, tool, ldt) result(rst)
    ! Constructs a mechanism and returns an opaque handle to it.
    logical, intent(in) :: planar
    integer(c_int), intent(in), value :: nlinks, njoints, base, effector, ldt
    type(c_mechanism_link), intent(in) :: links(nlinks)
    type(c_joint), intent(in) :: joints(njoints)
    type(c_ptr), intent(in), value :: tool
    type(c_ptr) :: rst

    integer(int32) :: i, nf
    real(real64) :: T(4, 4)
    real(real64), allocatable, dimension(:,:,:) :: frames
    real(c_double), pointer, dimension(:) :: fptr
    real(c_double), pointer, dimension(:,:) :: tptr
    type(link_container), allocatable, dimension(:) :: f_links
    type(joint), allocatable, dimension(:) :: f_joints
    type(c_mechanism_container), pointer :: cont
    type(parallel_linkage) :: spatial
    type(planar_linkage) :: flat

    ! Convert the links
    allocate(f_links(nlinks))
    do i = 1, nlinks
        nf = int(links(i)%frame_count, int32)
        if (nf < 1) error stop DYN_INVALID_INPUT_ERROR
        call c_f_pointer(links(i)%frames, fptr, [16 * nf])
        frames = reshape(fptr, [4, 4, nf])
        allocate(f_links(i)%item, source = multi_joint_link(frames, &
            links(i)%mass, reshape(links(i)%inertia, [3, 3]), links(i)%cg))
    end do

    ! Convert the joints
    allocate(f_joints(njoints))
    do i = 1, njoints
        f_joints(i) = joint( &
            int(joints(i)%joint_type, int32), &
            int(joints(i)%parent_link, int32), &
            int(joints(i)%child_link, int32), &
            int(joints(i)%parent_frame, int32), &
            int(joints(i)%child_frame, int32), &
            logical(joints(i)%actuated))
    end do

    ! The tool frame is optional; a null pointer denotes an identity matrix
    T = 0.0d0
    do i = 1, 4
        T(i,i) = 1.0d0
    end do
    if (c_associated(tool)) then
        if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
        call c_f_pointer(tool, tptr, [int(ldt, int32), 4])
        T = tptr(1:4,1:4)
    end if

    ! Construct the mechanism
    allocate(cont)
    if (planar) then
        flat = planar_linkage(f_links, f_joints, int(base, int32), &
            int(effector, int32), T)
        allocate(cont%item, source = flat)
    else
        spatial = parallel_linkage(f_links, f_joints, int(base, int32), &
            int(effector, int32), T)
        allocate(cont%item, source = spatial)
    end if
    rst = c_loc(cont)
end function

! ------------------------------------------------------------------------------
function c_create_parallel_linkage(nlinks, links, njoints, joints, base, &
    effector, tool, ldt) result(rst) &
    bind(C, name = "c_create_parallel_linkage")
    integer(c_int), intent(in), value :: nlinks, njoints, base, effector, ldt
    type(c_mechanism_link), intent(in) :: links(nlinks)
    type(c_joint), intent(in) :: joints(njoints)
    type(c_ptr), intent(in), value :: tool
    type(c_ptr) :: rst

    rst = build_mechanism(.false., nlinks, links, njoints, joints, base, &
        effector, tool, ldt)
end function

! ------------------------------------------------------------------------------
function c_create_planar_linkage(nlinks, links, njoints, joints, base, &
    effector, tool, ldt) result(rst) &
    bind(C, name = "c_create_planar_linkage")
    integer(c_int), intent(in), value :: nlinks, njoints, base, effector, ldt
    type(c_mechanism_link), intent(in) :: links(nlinks)
    type(c_joint), intent(in) :: joints(njoints)
    type(c_ptr), intent(in), value :: tool
    type(c_ptr) :: rst

    rst = build_mechanism(.true., nlinks, links, njoints, joints, base, &
        effector, tool, ldt)
end function

! ------------------------------------------------------------------------------
subroutine c_free_mechanism(obj) bind(C, name = "c_free_mechanism")
    type(c_ptr), intent(in), value :: obj

    type(c_mechanism_container), pointer :: cont

    if (.not.c_associated(obj)) return
    call c_f_pointer(obj, cont)
    if (allocated(cont%item)) deallocate(cont%item)
    deallocate(cont)
end subroutine

! ------------------------------------------------------------------------------
function c_mechanism_link_count(obj) result(rst) &
    bind(C, name = "c_mechanism_link_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_link_count(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_joint_count(obj) result(rst) &
    bind(C, name = "c_mechanism_joint_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_joint_count(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_variable_count(obj) result(rst) &
    bind(C, name = "c_mechanism_variable_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_variable_count(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_loop_count(obj) result(rst) &
    bind(C, name = "c_mechanism_loop_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_loop_count(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_constraint_count(obj) result(rst) &
    bind(C, name = "c_mechanism_constraint_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_constraint_count(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_degrees_of_freedom(obj) result(rst) &
    bind(C, name = "c_mechanism_degrees_of_freedom")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_degrees_of_freedom(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_actuated_variable_count(obj) result(rst) &
    bind(C, name = "c_mechanism_actuated_variable_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_actuated_variable_count(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_space_dimension(obj) result(rst) &
    bind(C, name = "c_mechanism_space_dimension")
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    mech => get_mechanism(obj)
    rst = 0
    if (associated(mech)) rst = int(mech%get_space_dimension(), c_int)
end function

! ------------------------------------------------------------------------------
function c_mechanism_link_frame_count(obj, i) result(rst) &
    bind(C, name = "c_mechanism_link_frame_count")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: i
    integer(c_int) :: rst

    class(kinematic_mechanism), pointer :: mech
    class(link), pointer :: lnk

    rst = 0
    mech => get_mechanism(obj)
    if (.not.associated(mech)) return
    lnk => mech%get_link(int(i, int32))
    rst = int(lnk%get_joint_count(), c_int)
end function

! ------------------------------------------------------------------------------
subroutine c_mechanism_link_frame(obj, i, k, T, ldt) &
    bind(C, name = "c_mechanism_link_frame")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: i, k, ldt
    real(c_double), intent(out) :: T(ldt, 4)

    class(kinematic_mechanism), pointer :: mech
    class(link), pointer :: lnk

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    lnk => mech%get_link(int(i, int32))
    T(1:4,1:4) = lnk%get_joint_frame(int(k, int32))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_get_configuration(obj, n, q) &
    bind(C, name = "c_mechanism_get_configuration")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n
    real(c_double), intent(out) :: q(n)

    class(kinematic_mechanism), pointer :: mech

    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    q = mech%get_configuration()
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_set_configuration(obj, n, q) &
    bind(C, name = "c_mechanism_set_configuration")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: q(n)

    class(kinematic_mechanism), pointer :: mech

    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    call mech%set_configuration(q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_body_transform(obj, i, n, q, T, ldt) &
    bind(C, name = "c_mechanism_body_transform")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: i, n, ldt
    real(c_double), intent(in) :: q(n)
    real(c_double), intent(out) :: T(ldt, 4)

    class(kinematic_mechanism), pointer :: mech

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    T(1:4,1:4) = mech%body_transform(int(i, int32), q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_end_effector_transform(obj, n, q, T, ldt) &
    bind(C, name = "c_mechanism_end_effector_transform")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n, ldt
    real(c_double), intent(in) :: q(n)
    real(c_double), intent(out) :: T(ldt, 4)

    class(kinematic_mechanism), pointer :: mech

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    T(1:4,1:4) = mech%end_effector_transform(q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_constraints(obj, n, q, nc, f) &
    bind(C, name = "c_mechanism_constraints")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n, nc
    real(c_double), intent(in) :: q(n)
    real(c_double), intent(out) :: f(nc)

    class(kinematic_mechanism), pointer :: mech

    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    f = mech%constraints(q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_constraint_jacobian(obj, n, q, jac, ldj) &
    bind(C, name = "c_mechanism_constraint_jacobian")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n, ldj
    real(c_double), intent(in) :: q(n)
    real(c_double), intent(out) :: jac(ldj, n)

    integer(int32) :: nc
    class(kinematic_mechanism), pointer :: mech

    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    nc = mech%get_constraint_count()
    if (ldj < nc) error stop DYN_INVALID_INPUT_ERROR
    jac(1:nc,1:n) = mech%constraint_jacobian(q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_solve_configuration(obj, na, qa, n, q, ib) &
    bind(C, name = "c_mechanism_solve_configuration")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: na, n
    real(c_double), intent(in) :: qa(na)
    real(c_double), intent(out) :: q(n)
    type(c_iteration_behavior), intent(out) :: ib

    class(kinematic_mechanism), pointer :: mech
    type(iteration_behavior) :: fib

    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    q = mech%solve_configuration(qa, ib = fib)
    ib = fib
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_forward_kinematics(obj, na, qa, T, ldt, ib) &
    bind(C, name = "c_mechanism_forward_kinematics")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: na, ldt
    real(c_double), intent(in) :: qa(na)
    real(c_double), intent(out) :: T(ldt, 4)
    type(c_iteration_behavior), intent(out) :: ib

    class(kinematic_mechanism), pointer :: mech
    type(iteration_behavior) :: fib

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    T(1:4,1:4) = mech%forward_kinematics(qa, ib = fib)
    ib = fib
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_jacobian(obj, na, qa, jac, ldj) &
    bind(C, name = "c_mechanism_jacobian")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: na, ldj
    real(c_double), intent(in) :: qa(na)
    real(c_double), intent(out) :: jac(ldj, na)

    integer(int32) :: nd
    class(kinematic_mechanism), pointer :: mech

    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    nd = mech%get_space_dimension()
    if (ldj < nd) error stop DYN_INVALID_INPUT_ERROR
    jac(1:nd,1:na) = mech%jacobian(qa)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_mechanism_inverse_kinematics(obj, trg, ldt, na, qa, ib) &
    bind(C, name = "c_mechanism_inverse_kinematics")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: ldt, na
    real(c_double), intent(in) :: trg(ldt, 4)
    real(c_double), intent(out) :: qa(na)
    type(c_iteration_behavior), intent(out) :: ib

    class(kinematic_mechanism), pointer :: mech
    type(iteration_behavior) :: fib

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    mech => get_mechanism(obj)
    if (.not.associated(mech)) error stop DYN_NULL_POINTER_ERROR
    qa = mech%inverse_kinematics(trg(1:4,1:4), ib = fib)
    ib = fib
end subroutine

end module
