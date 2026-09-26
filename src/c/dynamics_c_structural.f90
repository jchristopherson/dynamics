module dynamics_c_structural
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

pure function c_to_material(m) result(rst)
    ! Converts a c_material to a material.
    type(c_material), intent(in) :: m
    type(material) :: rst
    rst%density = m%density
    rst%modulus = m%modulus
    rst%poissons_ratio = m%poissons_ratio
end function

! ------------------------------------------------------------------------------
pure function c_to_node(n) result(rst)
    ! Converts a c_node to a node.
    type(c_node), intent(in) :: n
    type(node) :: rst
    rst%index = n%index
    rst%dof = n%dof
    rst%x = n%x
    rst%y = n%y
    rst%z = n%z
end function

! ------------------------------------------------------------------------------
pure function c_to_beam_element_2d(e) result(rst)
    ! Converts a c_beam_element_2d to a beam_element_2d.
    type(c_beam_element_2d), intent(in) :: e
    type(beam_element_2d) :: rst
    rst = beam_element_2d(c_to_material(e%material), e%area, &
        e%moment_of_inertia, c_to_node(e%node_1), c_to_node(e%node_2))
end function

! ------------------------------------------------------------------------------
pure function c_to_beam_element_3d(e) result(rst)
    ! Converts a c_beam_element_3d to a beam_element_3d.
    type(c_beam_element_3d), intent(in) :: e
    type(beam_element_3d) :: rst
    type(point) :: orient
    orient%x = e%orientation_point(1)
    orient%y = e%orientation_point(2)
    orient%z = e%orientation_point(3)
    rst = beam_element_3d(c_to_material(e%material), e%area, e%Ixx, e%Iyy, &
        e%Izz, e%Iyz, c_to_node(e%node_1), c_to_node(e%node_2), orient)
end function

! ------------------------------------------------------------------------------
function c_beam_element_2d_length(elem) result(rst) &
    bind(C, name = "c_beam_element_2d_length")
    type(c_beam_element_2d), intent(in) :: elem
    real(c_double) :: rst
    type(beam_element_2d) :: e
    e = c_to_beam_element_2d(elem)
    rst = e%length()
end function

! ------------------------------------------------------------------------------
subroutine c_beam_element_2d_stiffness_matrix(elem, rule, k, ldk) &
    bind(C, name = "c_beam_element_2d_stiffness_matrix")
    type(c_beam_element_2d), intent(in) :: elem
    integer(c_int), intent(in), value :: rule
    integer(c_int), intent(in), value :: ldk
    real(c_double), intent(out) :: k(ldk,6)
    type(beam_element_2d) :: e
    if (ldk < 6) error stop DYN_INVALID_INPUT_ERROR
    e = c_to_beam_element_2d(elem)
    if (rule > 0) then
        k(1:6,1:6) = e%stiffness_matrix(rule)
    else
        k(1:6,1:6) = e%stiffness_matrix()
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_2d_mass_matrix(elem, rule, m, ldm) &
    bind(C, name = "c_beam_element_2d_mass_matrix")
    type(c_beam_element_2d), intent(in) :: elem
    integer(c_int), intent(in), value :: rule
    integer(c_int), intent(in), value :: ldm
    real(c_double), intent(out) :: m(ldm,6)
    type(beam_element_2d) :: e
    if (ldm < 6) error stop DYN_INVALID_INPUT_ERROR
    e = c_to_beam_element_2d(elem)
    if (rule > 0) then
        m(1:6,1:6) = e%mass_matrix(rule)
    else
        m(1:6,1:6) = e%mass_matrix()
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_2d_rotation_matrix(elem, r, ldr) &
    bind(C, name = "c_beam_element_2d_rotation_matrix")
    type(c_beam_element_2d), intent(in) :: elem
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,6)
    type(beam_element_2d) :: e
    if (ldr < 6) error stop DYN_INVALID_INPUT_ERROR
    e = c_to_beam_element_2d(elem)
    r(1:6,1:6) = e%rotation_matrix()
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_2d_strain(elem, displacement, s, strain) &
    bind(C, name = "c_beam_element_2d_strain")
    type(c_beam_element_2d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(6)
    real(c_double), intent(in), value :: s
    real(c_double), intent(out) :: strain(2)
    type(beam_element_2d) :: e
    e = c_to_beam_element_2d(elem)
    strain = e%strain(displacement, [s])
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_2d_stress(elem, displacement, s, stress) &
    bind(C, name = "c_beam_element_2d_stress")
    type(c_beam_element_2d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(6)
    real(c_double), intent(in), value :: s
    real(c_double), intent(out) :: stress(2)
    type(beam_element_2d) :: e
    e = c_to_beam_element_2d(elem)
    stress = e%stress(displacement, [s])
end subroutine

! ------------------------------------------------------------------------------
function c_beam_element_2d_shear_force(elem, displacement, s) result(rst) &
    bind(C, name = "c_beam_element_2d_shear_force")
    type(c_beam_element_2d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(6)
    real(c_double), intent(in), value :: s
    real(c_double) :: rst
    type(beam_element_2d) :: e
    e = c_to_beam_element_2d(elem)
    rst = e%shear_force(displacement, [s])
end function

! ------------------------------------------------------------------------------
function c_beam_element_2d_bending_moment(elem, displacement, s) result(rst) &
    bind(C, name = "c_beam_element_2d_bending_moment")
    type(c_beam_element_2d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(6)
    real(c_double), intent(in), value :: s
    real(c_double) :: rst
    type(beam_element_2d) :: e
    e = c_to_beam_element_2d(elem)
    rst = e%bending_moment(displacement, [s])
end function

! ------------------------------------------------------------------------------
subroutine c_beam_element_2d_external_force_vector(elem, q, rule, f) &
    bind(C, name = "c_beam_element_2d_external_force_vector")
    type(c_beam_element_2d), intent(in) :: elem
    real(c_double), intent(in) :: q(2)
    integer(c_int), intent(in), value :: rule
    real(c_double), intent(out) :: f(6)
    type(beam_element_2d) :: e
    e = c_to_beam_element_2d(elem)
    if (rule > 0) then
        f = e%external_force_vector(q, rule)
    else
        f = e%external_force_vector(q)
    end if
end subroutine

! ------------------------------------------------------------------------------
function c_beam_element_3d_length(elem) result(rst) &
    bind(C, name = "c_beam_element_3d_length")
    type(c_beam_element_3d), intent(in) :: elem
    real(c_double) :: rst
    type(beam_element_3d) :: e
    e = c_to_beam_element_3d(elem)
    rst = e%length()
end function

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_stiffness_matrix(elem, rule, k, ldk) &
    bind(C, name = "c_beam_element_3d_stiffness_matrix")
    type(c_beam_element_3d), intent(in) :: elem
    integer(c_int), intent(in), value :: rule
    integer(c_int), intent(in), value :: ldk
    real(c_double), intent(out) :: k(ldk,12)
    type(beam_element_3d) :: e
    if (ldk < 12) error stop DYN_INVALID_INPUT_ERROR
    e = c_to_beam_element_3d(elem)
    if (rule > 0) then
        k(1:12,1:12) = e%stiffness_matrix(rule)
    else
        k(1:12,1:12) = e%stiffness_matrix()
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_mass_matrix(elem, rule, m, ldm) &
    bind(C, name = "c_beam_element_3d_mass_matrix")
    type(c_beam_element_3d), intent(in) :: elem
    integer(c_int), intent(in), value :: rule
    integer(c_int), intent(in), value :: ldm
    real(c_double), intent(out) :: m(ldm,12)
    type(beam_element_3d) :: e
    if (ldm < 12) error stop DYN_INVALID_INPUT_ERROR
    e = c_to_beam_element_3d(elem)
    if (rule > 0) then
        m(1:12,1:12) = e%mass_matrix(rule)
    else
        m(1:12,1:12) = e%mass_matrix()
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_rotation_matrix(elem, r, ldr) &
    bind(C, name = "c_beam_element_3d_rotation_matrix")
    type(c_beam_element_3d), intent(in) :: elem
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,12)
    type(beam_element_3d) :: e
    if (ldr < 12) error stop DYN_INVALID_INPUT_ERROR
    e = c_to_beam_element_3d(elem)
    r(1:12,1:12) = e%rotation_matrix()
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_strain(elem, displacement, s, strain) &
    bind(C, name = "c_beam_element_3d_strain")
    type(c_beam_element_3d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(12)
    real(c_double), intent(in), value :: s
    real(c_double), intent(out) :: strain(4)
    type(beam_element_3d) :: e
    e = c_to_beam_element_3d(elem)
    strain = e%strain(displacement, [s])
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_stress(elem, displacement, s, stress) &
    bind(C, name = "c_beam_element_3d_stress")
    type(c_beam_element_3d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(12)
    real(c_double), intent(in), value :: s
    real(c_double), intent(out) :: stress(4)
    type(beam_element_3d) :: e
    e = c_to_beam_element_3d(elem)
    stress = e%stress(displacement, [s])
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_shear_force(elem, displacement, s, force) &
    bind(C, name = "c_beam_element_3d_shear_force")
    type(c_beam_element_3d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(12)
    real(c_double), intent(in), value :: s
    real(c_double), intent(out) :: force(2)
    type(beam_element_3d) :: e
    e = c_to_beam_element_3d(elem)
    force = e%shear_force(displacement, [s])
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_bending_moment(elem, displacement, s, moment) &
    bind(C, name = "c_beam_element_3d_bending_moment")
    type(c_beam_element_3d), intent(in) :: elem
    real(c_double), intent(in) :: displacement(12)
    real(c_double), intent(in), value :: s
    real(c_double), intent(out) :: moment(3)
    type(beam_element_3d) :: e
    e = c_to_beam_element_3d(elem)
    moment = e%bending_moment(displacement, [s])
end subroutine

! ------------------------------------------------------------------------------
subroutine c_beam_element_3d_external_force_vector(elem, q, rule, f) &
    bind(C, name = "c_beam_element_3d_external_force_vector")
    type(c_beam_element_3d), intent(in) :: elem
    real(c_double), intent(in) :: q(4)
    integer(c_int), intent(in), value :: rule
    real(c_double), intent(out) :: f(12)
    type(beam_element_3d) :: e
    e = c_to_beam_element_3d(elem)
    if (rule > 0) then
        f = e%external_force_vector(q, rule)
    else
        f = e%external_force_vector(q)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine c_assemble_static_system_beam_2d(gdof, n, elements, nn, nodes, &
    rule, k, ldk) bind(C, name = "c_assemble_static_system_beam_2d")
    integer(c_int), intent(in), value :: gdof, n, nn, rule, ldk
    type(c_beam_element_2d), intent(in) :: elements(n)
    type(c_node), intent(in) :: nodes(nn)
    real(c_double), intent(out) :: k(ldk,gdof)

    integer(int32) :: i
    type(beam_element_2d), allocatable, dimension(:) :: felements
    type(node), allocatable, dimension(:) :: fnodes
    real(real64), allocatable, dimension(:,:) :: kf

    if (ldk < gdof) error stop DYN_INVALID_INPUT_ERROR

    allocate(felements(n), fnodes(nn))
    do i = 1, n
        felements(i) = c_to_beam_element_2d(elements(i))
    end do
    do i = 1, nn
        fnodes(i) = c_to_node(nodes(i))
    end do

    if (rule > 0) then
        call assemble_static_system(gdof, felements, fnodes, kf, rule)
    else
        call assemble_static_system(gdof, felements, fnodes, kf)
    end if
    k(1:gdof,1:gdof) = kf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_assemble_dynamic_system_beam_2d(gdof, n, elements, nn, nodes, &
    rule, m, ldm, k, ldk) bind(C, name = "c_assemble_dynamic_system_beam_2d")
    integer(c_int), intent(in), value :: gdof, n, nn, rule, ldm, ldk
    type(c_beam_element_2d), intent(in) :: elements(n)
    type(c_node), intent(in) :: nodes(nn)
    real(c_double), intent(out) :: m(ldm,gdof)
    real(c_double), intent(out) :: k(ldk,gdof)

    integer(int32) :: i
    type(beam_element_2d), allocatable, dimension(:) :: felements
    type(node), allocatable, dimension(:) :: fnodes
    real(real64), allocatable, dimension(:,:) :: mf, kf

    if (ldm < gdof) error stop DYN_INVALID_INPUT_ERROR
    if (ldk < gdof) error stop DYN_INVALID_INPUT_ERROR

    allocate(felements(n), fnodes(nn))
    do i = 1, n
        felements(i) = c_to_beam_element_2d(elements(i))
    end do
    do i = 1, nn
        fnodes(i) = c_to_node(nodes(i))
    end do

    if (rule > 0) then
        call assemble_dynamic_system(gdof, felements, fnodes, mf, kf, rule)
    else
        call assemble_dynamic_system(gdof, felements, fnodes, mf, kf)
    end if
    m(1:gdof,1:gdof) = mf
    k(1:gdof,1:gdof) = kf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_assemble_static_system_beam_3d(gdof, n, elements, nn, nodes, &
    rule, k, ldk) bind(C, name = "c_assemble_static_system_beam_3d")
    integer(c_int), intent(in), value :: gdof, n, nn, rule, ldk
    type(c_beam_element_3d), intent(in) :: elements(n)
    type(c_node), intent(in) :: nodes(nn)
    real(c_double), intent(out) :: k(ldk,gdof)

    integer(int32) :: i
    type(beam_element_3d), allocatable, dimension(:) :: felements
    type(node), allocatable, dimension(:) :: fnodes
    real(real64), allocatable, dimension(:,:) :: kf

    if (ldk < gdof) error stop DYN_INVALID_INPUT_ERROR

    allocate(felements(n), fnodes(nn))
    do i = 1, n
        felements(i) = c_to_beam_element_3d(elements(i))
    end do
    do i = 1, nn
        fnodes(i) = c_to_node(nodes(i))
    end do

    if (rule > 0) then
        call assemble_static_system(gdof, felements, fnodes, kf, rule)
    else
        call assemble_static_system(gdof, felements, fnodes, kf)
    end if
    k(1:gdof,1:gdof) = kf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_assemble_dynamic_system_beam_3d(gdof, n, elements, nn, nodes, &
    rule, m, ldm, k, ldk) bind(C, name = "c_assemble_dynamic_system_beam_3d")
    integer(c_int), intent(in), value :: gdof, n, nn, rule, ldm, ldk
    type(c_beam_element_3d), intent(in) :: elements(n)
    type(c_node), intent(in) :: nodes(nn)
    real(c_double), intent(out) :: m(ldm,gdof)
    real(c_double), intent(out) :: k(ldk,gdof)

    integer(int32) :: i
    type(beam_element_3d), allocatable, dimension(:) :: felements
    type(node), allocatable, dimension(:) :: fnodes
    real(real64), allocatable, dimension(:,:) :: mf, kf

    if (ldm < gdof) error stop DYN_INVALID_INPUT_ERROR
    if (ldk < gdof) error stop DYN_INVALID_INPUT_ERROR

    allocate(felements(n), fnodes(nn))
    do i = 1, n
        felements(i) = c_to_beam_element_3d(elements(i))
    end do
    do i = 1, nn
        fnodes(i) = c_to_node(nodes(i))
    end do

    if (rule > 0) then
        call assemble_dynamic_system(gdof, felements, fnodes, mf, kf, rule)
    else
        call assemble_dynamic_system(gdof, felements, fnodes, mf, kf)
    end if
    m(1:gdof,1:gdof) = mf
    k(1:gdof,1:gdof) = kf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_apply_boundary_conditions_mtx(n, nbc, gdof, x, ldx, rst, ldr) &
    bind(C, name = "c_apply_boundary_conditions_mtx")
    integer(c_int), intent(in), value :: n, nbc, ldx, ldr
    integer(c_int), intent(inout) :: gdof(nbc)
    real(c_double), intent(in) :: x(ldx,n)
    real(c_double), intent(out) :: rst(ldr, n - nbc)

    integer(int32), allocatable, dimension(:) :: fgdof
    real(real64), allocatable, dimension(:,:) :: frst

    if (ldx < n) error stop DYN_INVALID_INPUT_ERROR
    if (ldr < n - nbc) error stop DYN_INVALID_INPUT_ERROR

    fgdof = gdof
    frst = apply_boundary_conditions(fgdof, x(1:n,1:n))
    gdof = fgdof
    rst(1:n-nbc,1:n-nbc) = frst
end subroutine

! ------------------------------------------------------------------------------
subroutine c_apply_boundary_conditions_vec(n, nbc, gdof, x, rst) &
    bind(C, name = "c_apply_boundary_conditions_vec")
    integer(c_int), intent(in), value :: n, nbc
    integer(c_int), intent(inout) :: gdof(nbc)
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(out) :: rst(n - nbc)

    integer(int32), allocatable, dimension(:) :: fgdof
    real(real64), allocatable, dimension(:) :: frst

    fgdof = gdof
    frst = apply_boundary_conditions(fgdof, x)
    gdof = fgdof
    rst = frst
end subroutine

! ------------------------------------------------------------------------------
subroutine c_restore_constrained_values_dense(nred, nbc, gdof, x, rst) &
    bind(C, name = "c_restore_constrained_values_dense")
    integer(c_int), intent(in), value :: nred, nbc
    integer(c_int), intent(inout) :: gdof(nbc)
    real(c_double), intent(in) :: x(nred)
    real(c_double), intent(out) :: rst(nred + nbc)

    integer(int32), allocatable, dimension(:) :: fgdof
    real(real64), allocatable, dimension(:) :: frst

    fgdof = gdof
    frst = restore_constrained_values(fgdof, x)
    gdof = fgdof
    rst = frst
end subroutine

! ------------------------------------------------------------------------------
subroutine c_apply_displacement_constraint_dense(dof, val, n, k, ldk, f) &
    bind(C, name = "c_apply_displacement_constraint_dense")
    integer(c_int), intent(in), value :: dof, n, ldk
    real(c_double), intent(in), value :: val
    real(c_double), intent(inout) :: k(ldk,n)
    real(c_double), intent(inout) :: f(n)

    if (ldk < n) error stop DYN_INVALID_INPUT_ERROR

    call apply_displacement_constraint(dof, val, k(1:n,1:n), f)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_solve_static_system_dense(n, k, ldk, f, u) &
    bind(C, name = "c_solve_static_system_dense")
    integer(c_int), intent(in), value :: n, ldk
    real(c_double), intent(in) :: k(ldk,n)
    real(c_double), intent(in) :: f(n)
    real(c_double), intent(out) :: u(n)

    if (ldk < n) error stop DYN_INVALID_INPUT_ERROR

    u = solve_static_system(k(1:n,1:n), f)
end subroutine

! ------------------------------------------------------------------------------

end module
