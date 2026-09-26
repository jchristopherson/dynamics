module dynamics_c_linkage
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

! ------------------------------------------------------------------------------
subroutine c_build_serial_linkage(n, links, linkage) &
    bind(C, name = "c_build_serial_linkage")
    ! Arguments
    integer(c_int), intent(in), value :: n
    type(c_binary_link), intent(in) :: links(n)
    type(c_serial_linkage), intent(out) :: linkage

    ! Local Variables
    type(binary_link), allocatable, dimension(:) :: f_links
    type(serial_linkage) :: sf
    integer(int32) :: i, flag

    ! Process
    allocate(f_links(n))
    do i = 1, n
        f_links(i) = links(i)
    end do
    sf = serial_linkage(f_links)
    flag = c_alloc_serial_linkage(n, linkage)
    if (flag /= 0) return
    call convert_to_c_serial_linkage(linkage, sf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_serial_linkage_forward_kinematics(n, lnk, q, T, ldt) &
    bind(C, name = "c_serial_linkage_forward_kinematics")
    integer(c_int), intent(in), value :: n, ldt
    type(c_serial_linkage), intent(in) :: lnk
    real(c_double), intent(in) :: q(n)
    real(c_double), intent(out) :: T(ldt, 4)

    type(serial_linkage) :: f_lnk

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR

    f_lnk = lnk
    T(1:4,1:4) = f_lnk%forward_kinematics(q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_serial_linkage_jacobian(n, lnk, q, J, ldj) &
    bind(C, name = "c_serial_linkage_jacobian")
    integer(c_int), intent(in), value :: n, ldj
    type(c_serial_linkage), intent(in) :: lnk
    real(c_double), intent(in) :: q(n)
    real(c_double), intent(out) :: J(ldj, n)
    
    type(serial_linkage) :: f_lnk

    if (ldj < 6) error stop DYN_INVALID_INPUT_ERROR

    f_lnk = lnk
    J(1:6,1:n) = f_lnk%jacobian(q)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_serial_linkage_inverse_kinematics(n, lnk, qo, trg, ldt, q, ib) &
    bind(C, name = "c_serial_linkage_inverse_kinematics")
    integer(c_int), intent(in), value :: n, ldt
    type(c_serial_linkage), intent(in) :: lnk
    real(c_double), intent(in) :: qo(n)
    real(c_double), intent(in) :: trg(ldt, 4)
    real(c_double), intent(out) :: q(n)
    type(c_iteration_behavior), intent(out) :: ib

    type(serial_linkage) :: f_lnk
    type(iteration_behavior) :: fib

    if (ldt < 4) error stop DYN_INVALID_INPUT_ERROR

    f_lnk = lnk
    q = f_lnk%inverse_kinematics(qo, trg(1:4,1:4), fib)
    ib = fib
end subroutine

end module
