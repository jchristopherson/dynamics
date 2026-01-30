module dynamics_helper
    use iso_fortran_env
    implicit none
    private
    public :: cross_product
    public :: to_skew_symmetric
    public :: vector_angle
    public :: scalar_projection
    public :: vector_projection

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
pure function vector_angle(x, y) result(rst)
    !! Computes the angle between two vectors.
    real(real64), intent(in), dimension(:) :: x
        !! The first vector.
    real(real64), intent(in), dimension(size(x)) :: y
        !! The second vector.
    real(real64) :: rst
        !! The angle, in radians.

    ! Local Variables
    real(real64) :: xmag, ymag, ct

    ! Process
    xmag = norm2(x)
    ymag = norm2(y)
    ct = dot_product(x, y) / (xmag * ymag)
    rst = acos(ct)
end function

! ------------------------------------------------------------------------------
pure function scalar_projection(x, y) result(rst)
    !! Computes the projection of vector x onto vector y.  The scalar projection
    !! is defined such that \( s = \frac{\vec{x} \cdot \vec{y}}{||\vec{y}||} \).
    real(real64), intent(in), dimension(:) :: x
        !! The vector to project.
    real(real64), intent(in), dimension(size(x)) :: y
        !! The vector onto which x should be projected.
    real(real64) :: rst
        !! The scalar projection of x onto y.

    ! Process
    rst = dot_product(x, y) / norm2(y)
end function

! ------------------------------------------------------------------------------
pure function vector_projection(x, y) result(rst)
    !! Computes the vector pojection of vector x onto vector y.  The vector
    !! projection is defined such that \( proj_{y} \vec{x} =  
    !! \frac{\vec{x} \cdot \vec{y}}{||\vec{y}||^{2}} \vec{y}\)
    real(real64), intent(in), dimension(:) :: x
        !! The vector to project.
    real(real64), intent(in), dimension(size(x)) :: y
        !! The vector onto which x should be projected.
    real(real64) :: rst(3)
        !! The vector projection of x onto y.

    ! Process
    rst = y * dot_product(x, y) / dot_product(y, y)
end function

! ------------------------------------------------------------------------------
end module