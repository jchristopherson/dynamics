module dynamics_c_support
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

subroutine c_matmul(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) &
    bind(C, name = "c_matmul")
    use blas, only : DGEMM
    integer(c_int), intent(in), value :: m, n, k, lda, ldb, ldc
    real(c_double), intent(in), value :: alpha, beta
    real(c_double), intent(in) :: a(lda,k), b(ldb,n)
    real(c_double), intent(inout) :: c(ldc,n)
    call DGEMM('N', 'N', m, n, k, alpha, a(1:m,:), lda, b(1:k,:), ldb, beta, &
        c(1:m,:), ldc)
end subroutine

subroutine c_cross_product(x, y, z) bind(C, name = "c_cross_product")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double), intent(out) :: z(3)
    z = cross_product(x, y)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_to_skew_symmetric(x, y, ldy) bind(C, name = "c_to_skew_symmetric")
    real(c_double), intent(in) :: x(3)
    integer(c_int), intent(in), value :: ldy
    real(c_double), intent(out) :: y(ldy,3)
    if (ldy < 3) error stop DYN_INVALID_INPUT_ERROR
    y(1:3,1:3) = to_skew_symmetric(x)
end subroutine

! ------------------------------------------------------------------------------
function c_vector_angle(x, y) result(rst) bind(C, name = "c_vector_angle")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double) :: rst
    rst = vector_angle(x, y)
end function

! ------------------------------------------------------------------------------
function c_scalar_projection(x, y) result(rst) bind(C, name = "c_scalar_projection")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double) :: rst
    rst = scalar_projection(x, y)
end function

! ------------------------------------------------------------------------------
subroutine c_vector_projection(x, y, z) bind(C, name = "c_vector_projection")
    real(c_double), intent(in) :: x(3)
    real(c_double), intent(in) :: y(3)
    real(c_double), intent(out) :: z(3)
    z = vector_projection(x, y)
end subroutine

! ------------------------------------------------------------------------------
function c_vector_magnitude(n, x) result(rst) bind(C, name = "c_vector_magnitude")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double) :: rst
    rst = norm2(x)
end function

! ------------------------------------------------------------------------------
subroutine c_vector_normalize(n, x) bind(C, name = "c_vector_normalize")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(inout) :: x(n)
    x = x / norm2(x)
end subroutine

! ------------------------------------------------------------------------------
function c_dot_product(n, x, y) result(rst) bind(C, name = "c_dot_product")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double) :: rst
    rst = dot_product(x, y)
end function

end module
