module dynamics_c_variational
    use iso_c_binding
    use iso_fortran_env, only : int32, real64
    use dynamics
    use dynamics_error_handling
    use linalg, only : identity
    implicit none
    private

    type, bind(C) :: c_quaternion_vi
        real(c_double) :: w, x, y, z
    end type

    type, bind(C) :: c_rigid_body_vi
        real(c_double) :: mass
        real(c_double) :: cg(3)
        real(c_double) :: inertia(9)
    end type

    type, bind(C) :: c_variational_state_vi
        integer(c_int) :: body_count
        real(c_double) :: time
        type(c_ptr) :: position
        type(c_ptr) :: orientation
        type(c_ptr) :: velocity
        type(c_ptr) :: angular_velocity
    end type

    type, bind(C) :: c_variational_settings_vi
        real(c_double) :: tolerance
        real(c_double) :: finite_difference_step
        real(c_double) :: constraint_translation_scale
        real(c_double) :: constraint_rotation_scale
        integer(c_int) :: maximum_iterations
        integer(c_int) :: maximum_line_search_iterations
        integer(c_int) :: linear_solver
    end type

    type, bind(C) :: c_serial_linkage_vi
        integer(c_int) :: link_count
        type(c_ptr) :: links
    end type

    type, bind(C) :: c_binary_link_vi
        real(c_double) :: link_length, link_twist, link_offset, joint_angle
        integer(c_int) :: joint_type
        real(c_double) :: mass, cg(3), inertia(9)
    end type

    type, bind(C) :: c_mechanism_link_vi
        integer(c_int) :: frame_count
        type(c_ptr) :: frames
        real(c_double) :: mass, cg(3), inertia(9)
    end type

    type, bind(C) :: c_joint_vi
        integer(c_int) :: joint_type, parent_link, parent_frame
        integer(c_int) :: child_link, child_frame
        logical(c_bool) :: actuated
    end type

    type, bind(C) :: c_joint_reaction_vi
        real(c_double) :: force(3), moment(3)
    end type

    type, bind(C) :: c_linear_spring_vi
        integer(c_int) :: body_1, body_2
        real(c_double) :: point_1(3), point_2(3)
        real(c_double) :: stiffness, free_length
    end type

    type, bind(C) :: c_linear_damper_vi
        integer(c_int) :: body_1, body_2
        real(c_double) :: point_1(3), point_2(3)
        real(c_double) :: damping
    end type

    type, bind(C) :: c_torsional_spring_vi
        integer(c_int) :: joint_index
        real(c_double) :: stiffness, free_angle
    end type

    type, bind(C) :: c_torsional_damper_vi
        integer(c_int) :: joint_index
        real(c_double) :: damping
    end type

    type, bind(C) :: c_axial_element_result_vi
        real(c_double) :: length, length_rate, force
    end type

    type, bind(C) :: c_torsional_element_result_vi
        real(c_double) :: angle, angle_rate, torque
    end type

    abstract interface
        subroutine c_vi_force(state, force, torque, user_data) bind(C)
            import c_variational_state_vi, c_double, c_ptr
            type(c_variational_state_vi), intent(in) :: state
            real(c_double), intent(out) :: force(*), torque(*)
            type(c_ptr), intent(in), value :: user_data
        end subroutine
        subroutine c_vi_constraint(state, n, value, user_data) bind(C)
            import c_variational_state_vi, c_double, c_int, c_ptr
            type(c_variational_state_vi), intent(in) :: state
            integer(c_int), intent(in), value :: n
            real(c_double), intent(out) :: value(n)
            type(c_ptr), intent(in), value :: user_data
        end subroutine
        subroutine c_vi_jacobian(state, n, jac, ldj, user_data) bind(C)
            import c_variational_state_vi, c_double, c_int, c_ptr
            type(c_variational_state_vi), intent(in) :: state
            integer(c_int), intent(in), value :: n, ldj
            real(c_double), intent(out) :: jac(ldj,*)
            type(c_ptr), intent(in), value :: user_data
        end subroutine
        function c_vi_motion(t, user_data) result(rst) bind(C)
            import c_double, c_ptr
            real(c_double), intent(in), value :: t
            type(c_ptr), intent(in), value :: user_data
            real(c_double) :: rst
        end function
    end interface

    type c_callback_context
        procedure(c_vi_force), pointer, nopass :: force => null()
        procedure(c_vi_constraint), pointer, nopass :: constraint => null()
        procedure(c_vi_jacobian), pointer, nopass :: jacobian => null()
        type(c_ptr) :: user_data = c_null_ptr
    end type

    type c_dynamic_model_container
        type(linkage_dynamic_model), allocatable :: item
    end type

    procedure(c_vi_motion), pointer, save :: active_motion => null()
    type(c_ptr), save :: active_motion_data = c_null_ptr

    public :: c_vi_solve
    public :: c_vi_default_settings
    public :: c_vi_create_serial_model, c_vi_create_linkage_model
    public :: c_vi_free_linkage_model
    public :: c_vi_model_body_count, c_vi_model_joint_count
    public :: c_vi_model_constraint_count, c_vi_model_solve
    public :: c_vi_model_joint_reactions
    public :: c_vi_add_linear_spring, c_vi_add_linear_damper
    public :: c_vi_add_torsional_spring, c_vi_add_torsional_damper
    public :: c_vi_axial_count, c_vi_torsional_count
    public :: c_vi_axial_results, c_vi_torsional_results

contains

subroutine c_vi_default_settings(settings) &
    bind(C,name="c_default_variational_integrator_settings")
    type(c_variational_settings_vi), intent(out) :: settings
    type(variational_integrator_settings) :: defaults
    settings%tolerance = defaults%tolerance
    settings%finite_difference_step = defaults%finite_difference_step
    settings%constraint_translation_scale = defaults%constraint_translation_scale
    settings%constraint_rotation_scale = defaults%constraint_rotation_scale
    settings%maximum_iterations = defaults%maximum_iterations
    settings%maximum_line_search_iterations = defaults%maximum_line_search_iterations
    settings%linear_solver = defaults%linear_solver
end subroutine

pure function convert_settings(c) result(f)
    type(c_variational_settings_vi), intent(in) :: c
    type(variational_integrator_settings) :: f
    f%tolerance = c%tolerance
    f%finite_difference_step = c%finite_difference_step
    f%constraint_translation_scale = c%constraint_translation_scale
    f%constraint_rotation_scale = c%constraint_rotation_scale
    f%maximum_iterations = int(c%maximum_iterations, int32)
    f%maximum_line_search_iterations = int(c%maximum_line_search_iterations, int32)
    f%linear_solver = int(c%linear_solver, int32)
end function

subroutine pack_state(state, cstate, p, q, v, w)
    type(variational_state), intent(in) :: state
    type(c_variational_state_vi), intent(out) :: cstate
    real(c_double), allocatable, target, intent(out), dimension(:,:) :: p, v, w
    type(c_quaternion_vi), allocatable, target, intent(out), dimension(:) :: q
    integer(int32) :: i, n
    n = size(state%orientation)
    allocate(p(3,n), v(3,n), w(3,n), q(n))
    p = state%position; v = state%velocity; w = state%angular_velocity
    do i = 1, n
        q(i) = c_quaternion_vi(state%orientation(i)%w, state%orientation(i)%x, &
            state%orientation(i)%y, state%orientation(i)%z)
    end do
    cstate = c_variational_state_vi(int(n,c_int), state%time, c_loc(p), &
        c_loc(q), c_loc(v), c_loc(w))
end subroutine

subroutine force_bridge(t, state, force, torque, args)
    real(real64), intent(in) :: t
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:,:) :: force, torque
    class(*), intent(inout), optional :: args
    type(c_variational_state_vi) :: cs
    real(c_double), allocatable, target, dimension(:,:) :: p, v, w
    type(c_quaternion_vi), allocatable, target, dimension(:) :: q
    force = 0.0d0; torque = 0.0d0
    select type (ctx => args)
    type is (c_callback_context)
        if (.not.associated(ctx%force)) return
        call pack_state(state, cs, p, q, v, w)
        call ctx%force(cs, force, torque, ctx%user_data)
    end select
end subroutine

subroutine constraint_bridge(state, value, args)
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:) :: value
    class(*), intent(inout), optional :: args
    type(c_variational_state_vi) :: cs
    real(c_double), allocatable, target, dimension(:,:) :: p, v, w
    type(c_quaternion_vi), allocatable, target, dimension(:) :: q
    value = 0.0d0
    select type (ctx => args)
    type is (c_callback_context)
        if (.not.associated(ctx%constraint)) return
        call pack_state(state, cs, p, q, v, w)
        call ctx%constraint(cs, int(size(value),c_int), value, ctx%user_data)
    end select
end subroutine

subroutine jacobian_bridge(state, jac, args)
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:,:) :: jac
    class(*), intent(inout), optional :: args
    type(c_variational_state_vi) :: cs
    real(c_double), allocatable, target, dimension(:,:) :: p, v, w
    type(c_quaternion_vi), allocatable, target, dimension(:) :: q
    jac = 0.0d0
    select type (ctx => args)
    type is (c_callback_context)
        call pack_state(state, cs, p, q, v, w)
        call ctx%jacobian(cs, int(size(jac,1),c_int), jac, &
            int(size(jac,1),c_int), ctx%user_data)
    end select
end subroutine

subroutine copy_solution(solution, fm, p, q, v, w, m)
    type(variational_state), intent(in), dimension(:) :: solution
    real(real64), intent(in), dimension(:,:) :: fm
    real(c_double), intent(out), dimension(:,:,:) :: p, v, w
    type(c_quaternion_vi), intent(out), dimension(:,:) :: q
    real(c_double), intent(out), dimension(:,:) :: m
    integer(int32) :: i, j
    do j = 1, size(solution)
        p(:,:,j) = solution(j)%position
        v(:,:,j) = solution(j)%velocity
        w(:,:,j) = solution(j)%angular_velocity
        do i = 1, size(solution(j)%orientation)
            q(i,j) = c_quaternion_vi(solution(j)%orientation(i)%w, &
                solution(j)%orientation(i)%x, solution(j)%orientation(i)%y, &
                solution(j)%orientation(i)%z)
        end do
    end do
    if (size(m,1) > 0) m = fm
end subroutine

subroutine c_vi_solve(nbody, bodies, ntime, dt, p0, q0, v0, w0, nconstraint, &
    force_cb, constraint_cb, jacobian_cb, user_data, settings, p, q, v, w, m) &
    bind(C, name="c_variational_integrator_solve")
    integer(c_int), intent(in), value :: nbody, ntime, nconstraint
    real(c_double), intent(in), value :: dt
    type(c_rigid_body_vi), intent(in) :: bodies(nbody)
    real(c_double), intent(in) :: p0(3,nbody), v0(3,nbody), w0(3,nbody)
    type(c_quaternion_vi), intent(in) :: q0(nbody)
    type(c_funptr), intent(in), value :: force_cb, constraint_cb, jacobian_cb
    type(c_ptr), intent(in), value :: user_data
    type(c_variational_settings_vi), intent(in) :: settings
    real(c_double), intent(out) :: p(3,nbody,ntime), v(3,nbody,ntime), &
        w(3,nbody,ntime), m(nconstraint,ntime)
    type(c_quaternion_vi), intent(out) :: q(nbody,ntime)
    integer(int32) :: i
    type(rigid_body), allocatable, dimension(:) :: fb
    type(variational_state) :: state
    type(variational_state), allocatable, dimension(:) :: solution
    type(variational_integrator) :: integrator
    type(c_callback_context) :: ctx
    real(real64), allocatable, dimension(:,:) :: fm
    allocate(fb(nbody)); call initialize_variational_state(state, int(nbody,int32))
    do i = 1, nbody
        fb(i) = rigid_body(bodies(i)%mass, reshape(bodies(i)%inertia,[3,3]), &
            bodies(i)%cg)
        state%orientation(i) = quaternion([q0(i)%w,q0(i)%x,q0(i)%y,q0(i)%z])
    end do
    state%position=p0; state%velocity=v0; state%angular_velocity=w0
    integrator%settings = convert_settings(settings)
    if (c_associated(force_cb)) call c_f_procpointer(force_cb,ctx%force)
    if (c_associated(constraint_cb)) call c_f_procpointer(constraint_cb,ctx%constraint)
    if (c_associated(jacobian_cb)) call c_f_procpointer(jacobian_cb,ctx%jacobian)
    ctx%user_data=user_data
    if (associated(ctx%jacobian)) then
        solution=integrator%solve(fb,state,dt,int(ntime,int32), &
            constraint_count=int(nconstraint,int32),constraint=constraint_bridge, &
            force_function=force_bridge,constraint_jacobian=jacobian_bridge, &
            multipliers=fm,args=ctx)
    else
        solution=integrator%solve(fb,state,dt,int(ntime,int32), &
            constraint_count=int(nconstraint,int32),constraint=constraint_bridge, &
            force_function=force_bridge,multipliers=fm,args=ctx)
    end if
    call copy_solution(solution,fm,p,q,v,w,m)
end subroutine

function get_model(obj) result(model)
    type(c_ptr), intent(in), value :: obj
    type(linkage_dynamic_model), pointer :: model
    type(c_dynamic_model_container), pointer :: cont
    model=>null(); if (.not.c_associated(obj)) return
    call c_f_pointer(obj,cont); if (allocated(cont%item)) model=>cont%item
end function

function c_vi_create_serial_model(linkage,n,q) result(rst) &
    bind(C,name="c_create_serial_linkage_dynamic_model")
    type(c_serial_linkage_vi), intent(in) :: linkage
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: q(n)
    type(c_ptr) :: rst
    type(c_binary_link_vi), pointer, dimension(:) :: links
    type(binary_link), allocatable, dimension(:) :: flinks
    type(serial_linkage) :: serial
    type(c_dynamic_model_container), pointer :: cont
    integer(int32) :: i
    if (n /= linkage%link_count) error stop DYN_ARRAY_SIZE_ERROR
    call c_f_pointer(linkage%links,links,[linkage%link_count]); allocate(flinks(n))
    do i=1,n
        flinks(i)=binary_link(int(links(i)%joint_type,int32),links(i)%link_length, &
            links(i)%link_twist,links(i)%link_offset,links(i)%joint_angle, &
            links(i)%mass,reshape(links(i)%inertia,[3,3]),links(i)%cg)
    end do
    serial=serial_linkage(flinks); allocate(cont)
    allocate(cont%item,source=linkage_dynamic_model(serial,q)); rst=c_loc(cont)
end function

function c_vi_create_linkage_model(planar,nlinks,links,njoints,joints,base,nq,q) &
    result(rst) bind(C,name="c_create_linkage_dynamic_model")
    logical(c_bool), intent(in), value :: planar
    integer(c_int), intent(in), value :: nlinks,njoints,base,nq
    type(c_mechanism_link_vi), intent(in) :: links(nlinks)
    type(c_joint_vi), intent(in) :: joints(njoints)
    real(c_double), intent(in) :: q(nq)
    type(c_ptr) :: rst
    type(link_container), allocatable, dimension(:) :: flinks
    type(joint), allocatable, dimension(:) :: fjoints
    type(c_dynamic_model_container), pointer :: cont
    real(c_double), pointer, dimension(:) :: fp
    real(real64), allocatable, dimension(:,:,:) :: frames
    type(planar_linkage) :: pm
    type(parallel_linkage) :: sm
    integer(int32) :: i,nf
    allocate(flinks(nlinks),fjoints(njoints))
    do i=1,nlinks
        nf=links(i)%frame_count; call c_f_pointer(links(i)%frames,fp,[16*nf])
        frames=reshape(fp,[4,4,nf])
        allocate(flinks(i)%item,source=multi_joint_link(frames,links(i)%mass, &
            reshape(links(i)%inertia,[3,3]),links(i)%cg))
    end do
    do i=1,njoints
        fjoints(i)=joint(joints(i)%joint_type,joints(i)%parent_link, &
            joints(i)%child_link,joints(i)%parent_frame,joints(i)%child_frame, &
            logical(joints(i)%actuated))
    end do
    allocate(cont)
    if (planar) then
        pm=planar_linkage(flinks,fjoints,base=base)
        allocate(cont%item,source=linkage_dynamic_model(pm,q))
    else
        sm=parallel_linkage(flinks,fjoints,base=base)
        allocate(cont%item,source=linkage_dynamic_model(sm,q))
    end if
    rst=c_loc(cont)
end function

subroutine c_vi_free_linkage_model(obj) bind(C,name="c_free_linkage_dynamic_model")
    type(c_ptr), intent(in), value :: obj
    type(c_dynamic_model_container), pointer :: cont
    if (.not.c_associated(obj)) return; call c_f_pointer(obj,cont)
    if (allocated(cont%item)) deallocate(cont%item); deallocate(cont)
end subroutine

function c_vi_model_body_count(obj) result(rst) bind(C,name="c_linkage_dynamic_body_count")
    type(c_ptr), intent(in), value :: obj; integer(c_int) :: rst
    type(linkage_dynamic_model), pointer :: model
    model=>get_model(obj); rst=0; if(associated(model)) rst=model%get_body_count()
end function
function c_vi_model_joint_count(obj) result(rst) bind(C,name="c_linkage_dynamic_joint_count")
    type(c_ptr), intent(in), value :: obj; integer(c_int) :: rst
    type(linkage_dynamic_model), pointer :: model
    model=>get_model(obj); rst=0; if(associated(model)) rst=model%get_joint_count()
end function
function c_vi_model_constraint_count(obj) result(rst) bind(C,name="c_linkage_dynamic_constraint_count")
    type(c_ptr), intent(in), value :: obj; integer(c_int) :: rst
    type(linkage_dynamic_model), pointer :: model
    model=>get_model(obj); rst=0; if(associated(model)) rst=model%get_constraint_count()
end function

subroutine c_vi_add_linear_spring(obj,c) bind(C,name="c_linkage_dynamic_add_linear_spring")
    type(c_ptr), intent(in), value :: obj
    type(c_linear_spring_vi), intent(in) :: c
    type(linkage_dynamic_model), pointer :: model
    type(linear_spring) :: spring
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    spring%body_1=c%body_1; spring%body_2=c%body_2
    spring%point_1=c%point_1; spring%point_2=c%point_2
    spring%stiffness=c%stiffness; spring%free_length=c%free_length
    call model%add_linear_spring(spring)
end subroutine

subroutine c_vi_add_linear_damper(obj,c) bind(C,name="c_linkage_dynamic_add_linear_damper")
    type(c_ptr), intent(in), value :: obj
    type(c_linear_damper_vi), intent(in) :: c
    type(linkage_dynamic_model), pointer :: model
    type(linear_damper) :: damper
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    damper%body_1=c%body_1; damper%body_2=c%body_2
    damper%point_1=c%point_1; damper%point_2=c%point_2
    damper%damping=c%damping; call model%add_linear_damper(damper)
end subroutine

subroutine c_vi_add_torsional_spring(obj,c) bind(C,name="c_linkage_dynamic_add_torsional_spring")
    type(c_ptr), intent(in), value :: obj
    type(c_torsional_spring_vi), intent(in) :: c
    type(linkage_dynamic_model), pointer :: model
    type(torsional_spring) :: spring
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    spring%joint_index=c%joint_index; spring%stiffness=c%stiffness
    spring%free_angle=c%free_angle; call model%add_torsional_spring(spring)
end subroutine

subroutine c_vi_add_torsional_damper(obj,c) bind(C,name="c_linkage_dynamic_add_torsional_damper")
    type(c_ptr), intent(in), value :: obj
    type(c_torsional_damper_vi), intent(in) :: c
    type(linkage_dynamic_model), pointer :: model
    type(torsional_damper) :: damper
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    damper%joint_index=c%joint_index; damper%damping=c%damping
    call model%add_torsional_damper(damper)
end subroutine

function c_vi_axial_count(obj) result(rst) bind(C,name="c_linkage_dynamic_axial_element_count")
    type(c_ptr), intent(in), value :: obj; integer(c_int) :: rst
    type(linkage_dynamic_model), pointer :: model
    model=>get_model(obj); rst=0
    if(associated(model)) rst=model%get_axial_element_count()
end function

function c_vi_torsional_count(obj) result(rst) bind(C,name="c_linkage_dynamic_torsional_element_count")
    type(c_ptr), intent(in), value :: obj; integer(c_int) :: rst
    type(linkage_dynamic_model), pointer :: model
    model=>get_model(obj); rst=0
    if(associated(model)) rst=model%get_torsional_element_count()
end function

subroutine unpack_model_state(model,nbody,time,p,q,v,w,state)
    type(linkage_dynamic_model), intent(in) :: model
    integer(c_int), intent(in), value :: nbody
    real(c_double), intent(in), value :: time
    real(c_double), intent(in) :: p(3,nbody),v(3,nbody),w(3,nbody)
    type(c_quaternion_vi), intent(in) :: q(nbody)
    type(variational_state), intent(out) :: state
    integer(int32) :: i
    if(nbody/=model%get_body_count()) error stop DYN_ARRAY_SIZE_ERROR
    call initialize_variational_state(state,int(nbody,int32))
    state%time=time; state%position=p; state%velocity=v; state%angular_velocity=w
    do i=1,nbody
        state%orientation(i)=quaternion([q(i)%w,q(i)%x,q(i)%y,q(i)%z])
    end do
end subroutine

subroutine c_vi_axial_results(obj,nbody,time,p,q,v,w,r) &
    bind(C,name="c_linkage_dynamic_axial_element_results")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: nbody
    real(c_double), intent(in), value :: time
    real(c_double), intent(in) :: p(3,nbody),v(3,nbody),w(3,nbody)
    type(c_quaternion_vi), intent(in) :: q(nbody)
    type(c_axial_element_result_vi), intent(out) :: r(*)
    type(linkage_dynamic_model), pointer :: model
    type(variational_state) :: state
    type(axial_element_result), allocatable, dimension(:) :: fr
    integer(int32) :: i
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    call unpack_model_state(model,nbody,time,p,q,v,w,state)
    fr=model%get_axial_element_results(state)
    do i=1,size(fr)
        r(i)=c_axial_element_result_vi(fr(i)%length,fr(i)%length_rate,fr(i)%force)
    end do
end subroutine

subroutine c_vi_torsional_results(obj,nbody,time,p,q,v,w,r) &
    bind(C,name="c_linkage_dynamic_torsional_element_results")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: nbody
    real(c_double), intent(in), value :: time
    real(c_double), intent(in) :: p(3,nbody),v(3,nbody),w(3,nbody)
    type(c_quaternion_vi), intent(in) :: q(nbody)
    type(c_torsional_element_result_vi), intent(out) :: r(*)
    type(linkage_dynamic_model), pointer :: model
    type(variational_state) :: state
    type(torsional_element_result), allocatable, dimension(:) :: fr
    integer(int32) :: i
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    call unpack_model_state(model,nbody,time,p,q,v,w,state)
    fr=model%get_torsional_element_results(state)
    do i=1,size(fr)
        r(i)=c_torsional_element_result_vi(fr(i)%angle,fr(i)%angle_rate,fr(i)%torque)
    end do
end subroutine

function motion_bridge(t, args) result(rst)
    real(real64), intent(in) :: t
    class(*), intent(inout), optional :: args
    real(real64) :: rst
    rst=active_motion(t,active_motion_data)
end function

subroutine c_vi_model_solve(obj,nbody,nconstraint,settings,ntime,dt,gravity, &
    body_force,body_torque,prescribed_body,motion_cb,user_data,p,q,v,w,m) &
    bind(C,name="c_linkage_dynamic_solve")
    type(c_ptr), intent(in), value :: obj,user_data
    integer(c_int), intent(in), value :: nbody,nconstraint,ntime,prescribed_body
    real(c_double), intent(in), value :: dt
    type(c_variational_settings_vi), intent(in) :: settings
    real(c_double), intent(in) :: gravity(3),body_force(3,nbody),body_torque(3,nbody)
    type(c_funptr), intent(in), value :: motion_cb
    real(c_double), intent(out) :: p(3,nbody,ntime),v(3,nbody,ntime), &
        w(3,nbody,ntime),m(nconstraint,ntime)
    type(c_quaternion_vi), intent(out) :: q(nbody,ntime)
    type(linkage_dynamic_model), pointer :: model
    type(variational_integrator) :: integrator
    type(variational_state), allocatable, dimension(:) :: solution
    real(real64), allocatable, dimension(:,:) :: fm
    model=>get_model(obj); if(.not.associated(model)) error stop DYN_NULL_POINTER_ERROR
    if (nbody /= model%get_body_count()) error stop DYN_ARRAY_SIZE_ERROR
    if (nconstraint /= model%get_constraint_count() + &
        merge(1,0,c_associated(motion_cb))) error stop DYN_ARRAY_SIZE_ERROR
    integrator%settings=convert_settings(settings)
    if(c_associated(motion_cb)) then
        call c_f_procpointer(motion_cb,active_motion); active_motion_data=user_data
        solution=model%solve(integrator,dt,int(ntime,int32),gravity=gravity, &
            body_force=body_force,body_torque=body_torque, &
            prescribed_body=int(prescribed_body,int32), &
            prescribed_motion=motion_bridge,multipliers=fm)
        nullify(active_motion); active_motion_data=c_null_ptr
    else
        solution=model%solve(integrator,dt,int(ntime,int32),gravity=gravity, &
            body_force=body_force,body_torque=body_torque,multipliers=fm)
    end if
    if(size(fm,1)/=nconstraint) error stop DYN_ARRAY_SIZE_ERROR
    call copy_solution(solution,fm,p,q,v,w,m)
end subroutine

subroutine c_vi_model_joint_reactions(obj,nbody,nconstraint,time,p,q,v,w,m,r) &
    bind(C,name="c_linkage_dynamic_joint_reactions")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: nbody,nconstraint
    real(c_double), intent(in), value :: time
    real(c_double), intent(in) :: p(3,nbody),v(3,nbody),w(3,nbody),m(nconstraint)
    type(c_quaternion_vi), intent(in) :: q(nbody)
    type(c_joint_reaction_vi), intent(out) :: r(*)
    type(linkage_dynamic_model), pointer :: model
    type(variational_state) :: state
    type(joint_reaction), allocatable, dimension(:) :: fr
    integer(int32) :: i
    model=>get_model(obj); call initialize_variational_state(state,nbody)
    state%time=time; state%position=p; state%velocity=v; state%angular_velocity=w
    do i=1,nbody
        state%orientation(i)=quaternion([q(i)%w,q(i)%x,q(i)%y,q(i)%z])
    end do
    fr=model%get_joint_reactions(state,m)
    do i=1,size(fr)
        r(i)%force=fr(i)%force; r(i)%moment=fr(i)%moment
    end do
end subroutine

end module
