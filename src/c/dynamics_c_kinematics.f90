module dynamics_c_kinematics
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

subroutine c_dh_forward_kinematics_table(tbl, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_table")
    type(c_dh_table), intent(in) :: tbl
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt, 4)
    type(dh_table) :: ftbl
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    ftbl = tbl
    T(1:4,1:4) = dh_forward_kinematics(ftbl)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics(n, alpha, a, theta, d, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: alpha(n)
    real(c_double), intent(in) :: a(n)
    real(c_double), intent(in) :: theta(n)
    real(c_double), intent(in) :: d(n)
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(alpha, a, theta, d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_2(T1, ldt1, T2, ldt2, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_2")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_3(T1, ldt1, T2, ldt2, T3, ldt3, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_3")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt3 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_4(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_4")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt3 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt4 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_5(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_5")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt3 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt4 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt5 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_6(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T6, ldt6, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_6")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt6
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(in) :: T6(ldt6,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt3 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt4 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt5 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt6 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4), T6(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_7(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T6, ldt6, T7, ldt7, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_7")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt6
    integer(c_int), intent(in), value :: ldt7
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(in) :: T6(ldt6,4)
    real(c_double), intent(in) :: T7(ldt7,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt3 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt4 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt5 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt6 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt7 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4), T6(1:4,1:4), T7(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_forward_kinematics_8(T1, ldt1, T2, ldt2, T3, ldt3, T4, ldt4, &
    T5, ldt5, T6, ldt6, T7, ldt7, T8, ldt8, T, ldt) &
    bind(C, name = "c_dh_forward_kinematics_8")
    integer(c_int), intent(in), value :: ldt1
    integer(c_int), intent(in), value :: ldt2
    integer(c_int), intent(in), value :: ldt3
    integer(c_int), intent(in), value :: ldt4
    integer(c_int), intent(in), value :: ldt5
    integer(c_int), intent(in), value :: ldt6
    integer(c_int), intent(in), value :: ldt7
    integer(c_int), intent(in), value :: ldt8
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(in) :: T1(ldt1,4)
    real(c_double), intent(in) :: T2(ldt2,4)
    real(c_double), intent(in) :: T3(ldt3,4)
    real(c_double), intent(in) :: T4(ldt4,4)
    real(c_double), intent(in) :: T5(ldt5,4)
    real(c_double), intent(in) :: T6(ldt6,4)
    real(c_double), intent(in) :: T7(ldt7,4)
    real(c_double), intent(in) :: T8(ldt8,4)
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt1 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt2 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt3 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt4 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt5 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt6 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt7 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt8 < 4) error stop DYN_INVALID_INPUT_ERROR
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_forward_kinematics(T1(1:4,1:4), T2(1:4,1:4), T3(1:4,1:4), &
        T4(1:4,1:4), T5(1:4,1:4), T6(1:4,1:4), T7(1:4,1:4), T8(1:4,1:4))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_jacobian(n, alpha, a, theta, d, jtypes, jac, ldjac) &
    bind(C, name = "c_dh_jacobian")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: ldjac
    real(c_double), intent(in) :: alpha(n)
    real(c_double), intent(in) :: a(n)
    real(c_double), intent(in) :: theta(n)
    real(c_double), intent(in) :: d(n)
    integer(c_int), intent(in) :: jtypes(n)
    real(c_double), intent(out) :: jac(ldjac,n)
    if (ldjac < 6) error stop DYN_INVALID_INPUT_ERROR
    jac(1:6,1:n) = dh_jacobian(alpha, a, theta, d, jtypes)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_matrix(alpha, a, theta, d, T, ldt) &
    bind(C, name = "c_dh_matrix")
    real(c_double), intent(in), value :: alpha
    real(c_double), intent(in), value :: a
    real(c_double), intent(in), value :: theta
    real(c_double), intent(in), value :: d
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_matrix(alpha, a, theta, d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_rotate_x(alpha, T, ldt) &
    bind(C, name = "c_dh_rotate_x")
    real(c_double), intent(in), value :: alpha
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_rotate_x(alpha)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_rotate_z(theta, T, ldt) &
    bind(C, name = "c_dh_rotate_z")
    real(c_double), intent(in), value :: theta
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_rotate_z(theta)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_translate_x(a, T, ldt) &
    bind(C, name = "c_dh_translate_x")
    real(c_double), intent(in), value :: a
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_translate_x(a)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_dh_translate_z(d, T, ldt) &
    bind(C, name = "c_dh_translate_z")
    real(c_double), intent(in), value :: d
    integer(c_int), intent(in), value :: ldt
    real(c_double), intent(out) :: T(ldt,4)
    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR
    T(1:4,1:4) = dh_translate_z(d)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_jacobian_generating_vector(d, k, R, ldr, jtype, jvec) &
    bind(C, name = "c_jacobian_generating_vector")
    real(c_double), intent(in) :: d(3)
    real(c_double), intent(in) :: k(3)
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(in) :: R(ldr,3)
    integer(c_int), intent(in), value :: jtype
    real(c_double), intent(out) :: jvec(6)
    if (ldr < 3) error stop DYN_INVALID_INPUT_ERROR
    jvec = jacobian_generating_vector(d, k, R(1:3,1:3), jtype)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_solve_inverse_kinematics(njoints, neqn, mdl, qo, constraints, &
    qmax, qmin, jvar, resid, ib) &
    bind(C, name = "c_solve_inverse_kinematics")
    integer(c_int), intent(in), value :: njoints
    integer(c_int), intent(in), value :: neqn
    type(c_funptr), intent(in), value :: mdl
    real(c_double), intent(in) :: qo(njoints)
    real(c_double), intent(in) :: constraints(neqn)
    real(c_double), intent(in) :: qmax(njoints)
    real(c_double), intent(in) :: qmin(njoints)
    real(c_double), intent(out) :: jvar(njoints)
    real(c_double), intent(out) :: resid(neqn)
    type(c_iteration_behavior), intent(out) :: ib
    type(iteration_behavior) :: iter
    procedure(vecfcn), pointer :: fcn
    procedure(c_vecfcn), pointer :: fptr
    type(c_vecfcn_container) :: arg
    call c_f_procpointer(mdl, fptr)
    fcn => sik_fcn
    arg%fcn => fptr
    
    jvar = solve_inverse_kinematics(fcn, qo, constraints, df = resid, &
        qmax = qmax, qmin = qmin, ib = iter, args = arg)
    ib = iter
end subroutine

! --------------------
subroutine sik_fcn(x, f, args)
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args
    select type (args)
    class is (c_vecfcn_container)
        call args%fcn(size(x), size(f), x, f)
    end select
end subroutine


subroutine c_define_link_csys(xim1, zim1, zi, rim1, ri, csys) &
    bind(C, name = "c_define_link_csys")
    real(c_double), intent(in) :: xim1(3)
    real(c_double), intent(in) :: zim1(3)
    real(c_double), intent(in) :: zi(3)
    real(c_double), intent(in) :: rim1(3)
    real(c_double), intent(in) :: ri(3)
    type(c_coordinate_system), intent(out) :: csys
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_define_csys(i, j, k, o, csys) bind(C, name = "c_define_csys")
    real(c_double), intent(in) :: i(3)
    real(c_double), intent(in) :: j(3)
    real(c_double), intent(in) :: k(3)
    real(c_double), intent(in) :: o(3)
    type(c_coordinate_system), intent(out) :: csys
    csys = coordinate_system(i, j, k, o)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_build_dh_table(n, csys, tbl) bind(C, name = "c_build_dh_table")
    integer(c_int), intent(in), value :: n
    type(c_coordinate_system), intent(in) :: csys(n)
    type(c_dh_table), intent(out) :: tbl
    type(dh_table) :: ftbl
    integer(int32) :: i, flag
    type(coordinate_system), allocatable, dimension(:) :: c
    allocate(c(n))
    do i = 1, n
        c(i) = csys(i)
    end do
    ftbl = dh_table(c)
    flag = c_alloc_dh_table(size(ftbl%parameters), tbl)
    if (flag /= 0) return
    call convert_to_c_dh_table(tbl, ftbl)
end subroutine

end module
