module dynamics_helper
    use iso_fortran_env
    implicit none
    private
    public :: cross_product

contains
! ------------------------------------------------------------------------------
pure function cross_product(x, y) result(rst)
    !! Computes the cross-product of a vector.
    real(real64), intent(in) :: x(3)
        !! The left-hand-side argument.
    real(real64), intent(in) :: y(3)
        !! The right-hand-side argument
    real(real64) :: rst(3)
        !! The resulting vector.

    rst(1) = x(2) * y(3) - x(3) * y(2)
    rst(2) = x(3) * y(1) - x(1) * y(3)
    rst(3) = x(1) * y(2) - x(2) * y(1)
end function

! ------------------------------------------------------------------------------
end module