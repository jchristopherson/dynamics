module dynamics_helper
    use iso_fortran_env
    implicit none
    private
    public :: cross_product
    public :: to_skew_symmetric

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
pure function to_skew_symmetric(x) result(rst)
    !! Converts a 3-element vector to a 3-by-3 skew-symmetric matrix.  A 
    !! skew-symmetric matrix is defined as follows.
    !!
    !! $$ \tilde{x} = \left[ \begin{matrix} 0 & -x_{3} & x_{2} \\
    !! x_{3} & 0 & -x_{1} \\ -x_{2} & x_{1} & 0 \end{matrix} \right] $$
    real(real64), intent(in) :: x(3)
        !! The vector.
    real(real64) :: rst(3, 3)
        !! The resulting skew-symmetric matrix.

    ! Process
    rst = reshape([ &
        0.0d0, x(3), -x(2), &
        -x(3), 0.0d0, x(1), &
        x(2), -x(1), 0.0d0 &
    ], [3, 3])
end function

! ------------------------------------------------------------------------------
end module