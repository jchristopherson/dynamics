module dynamics_c_stability
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

subroutine c_determine_local_stability(n, a, lda, ev, flag) &
    bind(C, name = "c_determine_local_stability")
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: lda
    real(c_double), intent(in) :: a(lda,n)
    complex(c_double), intent(out) :: ev(n)
    integer(c_int), intent(out) :: flag
    
    if (lda < n) error stop DYN_INVALID_INPUT_ERROR
    flag = determine_local_stability(a(1:n,1:n), ev = ev)
end subroutine

end module
