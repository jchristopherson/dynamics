! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
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
    !! Computes the cross-product of two three-dimensional vectors.
    !! The result is orthogonal to both inputs and is defined by
    !! $$ \boldsymbol{x}\times\boldsymbol{y} =
    !! \begin{bmatrix}x_2y_3-x_3y_2\\x_3y_1-x_1y_3\\x_1y_2-x_2y_1\end{bmatrix}. $$
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
    !! Computes the unsigned angle between two nonzero vectors.
    !! $$ \theta = \cos^{-1}\left(\frac{\boldsymbol{x}\cdot\boldsymbol{y}}
    !! {\|\boldsymbol{x}\|\,\|\boldsymbol{y}\|}\right),\qquad 0\leq\theta\leq\pi. $$
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
    !! is defined such that
    !! $$ s = \frac{\boldsymbol{x}\cdot\boldsymbol{y}}{\|\boldsymbol{y}\|}. $$
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
    !! Computes the vector projection of vector x onto vector y.  The vector
    !! projection is defined such that \( proj_{y} \vec{x} =  
    !! $$ \operatorname{proj}_{\boldsymbol{y}}\boldsymbol{x} =
    !! \frac{\boldsymbol{x}\cdot\boldsymbol{y}}{\|\boldsymbol{y}\|^{2}}
    !! \boldsymbol{y}. $$
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