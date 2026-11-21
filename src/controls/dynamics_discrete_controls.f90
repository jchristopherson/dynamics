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

!> Defines discrete-time state-space and transfer function models along with
!! conversions between continuous and discrete forms.
module dynamics_discrete_controls
    use iso_fortran_env
    use nonlin_polynomials
    use linalg, only : mtx_inverse, eigen, identity
    use dynamics_error_handling, only : DYN_INVALID_INPUT_ERROR
    use dynamics_helper, only : matrix_exponential, siso_resolvent_coefficients
    use dynamics_controls, only : state_space, transfer_function
    implicit none
    private
    public :: discrete_state_space
    public :: discrete_transfer_function
    public :: c2d
    public :: d2c
    public :: DISCRETIZE_ZOH
    public :: DISCRETIZE_TUSTIN

    integer, parameter :: DISCRETIZE_ZOH = 1
        !! Zero-order hold discretization.
    integer, parameter :: DISCRETIZE_TUSTIN = 2
        !! Bilinear (Tustin) discretization.
    integer, parameter :: MAX_SQRT_ITERATIONS = 100
    integer, parameter :: MAX_SERIES_TERMS = 60

! ------------------------------------------------------------------------------
    type discrete_state_space
        !! Defines a discrete-time state-space model of the form
        !! $$ x_{k+1}=Ax_k+Bu_k,\qquad y_k=Cx_k+Du_k. $$
        real(real64), allocatable, dimension(:,:) :: A
            !! The N-by-N state matrix.
        real(real64), allocatable, dimension(:,:) :: B
            !! The N-by-M input matrix.
        real(real64), allocatable, dimension(:,:) :: C
            !! The P-by-N output matrix.
        real(real64), allocatable, dimension(:,:) :: D
            !! The P-by-M feedthrough matrix.
        real(real64) :: sample_time
            !! The sample time, \(T_s > 0\), of the model.
    contains
        procedure, public :: poles => dss_poles
        procedure, public :: to_transfer_function => dss_to_tf
    end type

    interface discrete_state_space
        module procedure :: dss_init
    end interface

! ------------------------------------------------------------------------------
    type discrete_transfer_function
        !! Defines a discrete-time transfer function of the form
        !! $$ H(z) = \frac{Y(z^{-1})}{X(z^{-1})}. $$
        !! The polynomial coefficients are stored in ascending powers of
        !! \(z^{-1}\), such that \( y_1 + y_2 z^{-1} + y_3 z^{-2} ... \).
        type(polynomial) :: Y
            !! The numerator polynomial in \(z^{-1}\).
        type(polynomial) :: X
            !! The denominator polynomial in \(z^{-1}\).
        real(real64) :: sample_time
            !! The sample time, \(T_s > 0\), of the model.
    contains
        procedure, public :: poles => dtf_poles
        procedure, public :: zeros => dtf_zeros
        procedure, public :: evaluate => dtf_eval_omega
        procedure, public :: to_ccf_state_space => dtf_to_ccf_statespace
        procedure, public :: to_ocf_state_space => dtf_to_ocf_statespace
    end type

    interface discrete_transfer_function
        module procedure :: dtf_init
    end interface

contains
! ------------------------------------------------------------------------------
pure function dss_init(a, b, c, d, sample_time) result(rst)
    !! Initializes a discrete-time state-space model.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N state matrix.
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-M input matrix.
    real(real64), intent(in), dimension(:,:) :: c
        !! The P-by-N output matrix.
    real(real64), intent(in), dimension(:,:) :: d
        !! The P-by-M feedthrough matrix.
    real(real64), intent(in) :: sample_time
        !! The sample time, which must be positive.
    type(discrete_state_space) :: rst
        !! The resulting [[discrete_state_space]].

    ! Local Variables
    type(state_space) :: ss

    ! Process
    if (sample_time <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    ss = state_space(a, b, c, d)
    rst%A = ss%A
    rst%B = ss%B
    rst%C = ss%C
    rst%D = ss%D
    rst%sample_time = sample_time
end function

! ------------------------------------------------------------------------------
pure function dtf_init(y, x, sample_time) result(rst)
    !! Initializes a discrete-time transfer function.
    real(real64), intent(in), dimension(:) :: y
        !! The numerator coefficients in ascending powers of \(z^{-1}\).
    real(real64), intent(in), dimension(:) :: x
        !! The denominator coefficients in ascending powers of \(z^{-1}\).
    real(real64), intent(in) :: sample_time
        !! The sample time, which must be positive.
    type(discrete_transfer_function) :: rst
        !! The resulting [[discrete_transfer_function]].

    ! Process
    if (size(y) < 1 .or. size(x) < 1) error stop DYN_INVALID_INPUT_ERROR
    if (sample_time <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    call rst%Y%initialize(y)
    call rst%X%initialize(x)
    rst%sample_time = sample_time
end function

! ------------------------------------------------------------------------------
pure function dss_poles(this) result(rst)
    !! Computes the poles of the discrete-time state-space model.  The poles
    !! are the eigenvalues of \(A\), and are the roots of \(\det(zI-A)=0\).
    class(discrete_state_space), intent(in) :: this
        !! The [[discrete_state_space]] object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The poles of the model.

    ! Local Variables
    integer :: n

    ! Process
    n = size(this%A, 1)
    allocate(rst(n))
    if (n == 0) return
    call eigen(this%A, rst)
end function

! ------------------------------------------------------------------------------
pure function dss_to_tf(this) result(rst)
    !! Converts a single-input, single-output discrete-time state-space model
    !! into a transfer function.  The result is expressed in ascending powers
    !! of \(z^{-1}\) by multiplying the numerator and denominator of
    !! $$ H(z)=C(zI-A)^{-1}B+D $$
    !! by \(z^{-n}\).
    class(discrete_state_space), intent(in) :: this
        !! The [[discrete_state_space]] object.
    type(discrete_transfer_function) :: rst
        !! The resulting [[discrete_transfer_function]].

    ! Local Variables
    integer :: i, n
    real(real64) :: d
    real(real64), allocatable, dimension(:) :: den, markov, num

    ! Process
    n = size(this%A, 1)
    d = this%D(1,1)
    call siso_resolvent_coefficients(this%A, this%B, this%C, den, markov)
    allocate(num(n + 1))
    num(1) = d * den(1)
    do i = 2, n + 1
        num(i) = d * den(i) + markov(i - 1)
    end do
    rst = discrete_transfer_function(num, den, this%sample_time)
end function

! ------------------------------------------------------------------------------
pure function dtf_poles(this) result(rst)
    !! Computes the poles of the discrete-time transfer function.
    class(discrete_transfer_function), intent(in) :: this
        !! The [[discrete_transfer_function]] object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The poles of the transfer function.

    ! Process
    rst = z_roots(this%X%get_all())
end function

! ------------------------------------------------------------------------------
pure function dtf_zeros(this) result(rst)
    !! Computes the zeros of the discrete-time transfer function.  Zeros at
    !! \(z=\infty\) that arise from a numerator with fewer coefficients than
    !! the denominator are not returned.
    class(discrete_transfer_function), intent(in) :: this
        !! The [[discrete_transfer_function]] object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The zeros of the transfer function.

    ! Process
    rst = z_roots(this%Y%get_all())
end function

! ------------------------------------------------------------------------------
pure elemental function dtf_eval_omega(this, omega) result(rst)
    !! Evaluates the transfer function at the specified frequency.
    class(discrete_transfer_function), intent(in) :: this
        !! The [[discrete_transfer_function]] object.
    real(real64), intent(in) :: omega
        !! The frequency, in rad/s.
    complex(real64) :: rst
        !! The result.

    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    integer :: j, nb, na
    complex(real64) :: delay, num, den

    na = this%X%order() + 1
    nb = this%Y%order() + 1
    delay = exp(cmplx(0.0d0, -omega * this%sample_time, kind = real64))
    num = zero
    do j = 1, nb
        num = num + this%Y%get(j) * delay**(j - 1)
    end do
    den = zero
    do j = 1, na
        den = den + this%X%get(j) * delay**(j - 1)
    end do
    rst = num / den
end function

! ------------------------------------------------------------------------------
pure function dtf_to_ccf_statespace(this) result(rst)
    !! Converts a causal discrete transfer function to controllable canonical
    !! state-space form.  The denominator's constant term in \(z^{-1}\)
    !! must be nonzero.
    class(discrete_transfer_function), intent(in) :: this
        !! The transfer function to convert.
    type(discrete_state_space) :: rst
        !! The resulting state-space model, with the original sample time.

    rst = dtf_to_canonical(this, .false.)
end function

! ------------------------------------------------------------------------------
pure function dtf_to_ocf_statespace(this) result(rst)
    !! Converts a causal discrete transfer function to observable canonical
    !! state-space form.  The denominator's constant term in \(z^{-1}\)
    !! must be nonzero.
    class(discrete_transfer_function), intent(in) :: this
        !! The transfer function to convert.
    type(discrete_state_space) :: rst
        !! The resulting state-space model, with the original sample time.

    rst = dtf_to_canonical(this, .true.)
end function

! ------------------------------------------------------------------------------
pure function dtf_to_canonical(this, observable) result(rst)
    class(discrete_transfer_function), intent(in) :: this
    logical, intent(in) :: observable
    type(discrete_state_space) :: rst

    integer :: i
    integer :: na
    integer :: nb
    integer :: n
    real(real64), allocatable, dimension(:) :: a
    real(real64), allocatable, dimension(:) :: b
    type(transfer_function) :: tf
    type(state_space) :: ss

    na = this%X%order() + 1
    nb = this%Y%order() + 1
    if (na < 1 .or. nb < 1) error stop DYN_INVALID_INPUT_ERROR
    if (this%sample_time <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    if (this%X%get(1) == 0.0d0) error stop DYN_INVALID_INPUT_ERROR

    ! Pad in z**(-1), then reverse to obtain polynomials in z.  A numerator
    ! longer than the denominator represents delays, not an improper model.
    do while (na > 1)
        if (this%X%get(na) /= 0.0d0) exit
        na = na - 1
    end do
    do while (nb > 1)
        if (this%Y%get(nb) /= 0.0d0) exit
        nb = nb - 1
    end do
    n = max(na, nb) - 1
    allocate(a(n + 1), b(n + 1), source = 0.0d0)
    do i = 1, na
        a(i) = this%X%get(i)
    end do
    do i = 1, nb
        b(i) = this%Y%get(i)
    end do
    tf = transfer_function(b(n+1:1:-1), a(n+1:1:-1))
    if (observable) then
        ss = tf%to_ocf_state_space()
    else
        ss = tf%to_ccf_state_space()
    end if
    rst%A = ss%A
    rst%B = ss%B
    rst%C = ss%C
    rst%D = ss%D
    rst%sample_time = this%sample_time
end function

! ------------------------------------------------------------------------------
pure function c2d(sys, sample_time, method) result(rst)
    !! Converts a continuous-time state-space model into a discrete-time
    !! state-space model using either a zero-order hold (ZOH) or the bilinear
    !! (Tustin) transformation.
    !! The ZOH method uses the exponential of the augmented matrix
    !! $$ \exp\left(\begin{bmatrix}A&B\\0&0\end{bmatrix}T_s\right) $$
    !! to obtain \(A_d\) and \(B_d\).  The Tustin method is
    !! $$ A_d=F(I+\tfrac{T_s}{2}A),\quad B_d=T_sF^2B,\quad
    !!    D_d=D+\tfrac{T_s}{2}CFB,\quad F=(I-\tfrac{T_s}{2}A)^{-1}. $$
    type(state_space), intent(in) :: sys
        !! The continuous-time [[state_space]] model.
    real(real64), intent(in) :: sample_time
        !! The sample time, which must be positive.
    integer, intent(in), optional :: method
        !! The discretization method, either [[DISCRETIZE_ZOH]] (the default)
        !! or [[DISCRETIZE_TUSTIN]].
    type(discrete_state_space) :: rst
        !! The resulting [[discrete_state_space]].

    ! Local Variables
    integer :: m, method_, n
    real(real64), allocatable, dimension(:,:) :: ad, bd, dd, f, fb, ident, mtx, phi

    ! Process
    if (sample_time <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    method_ = DISCRETIZE_ZOH
    if (present(method)) method_ = method
    n = size(sys%A, 1)
    m = size(sys%B, 2)
    select case (method_)
    case (DISCRETIZE_ZOH)
        allocate(mtx(n + m, n + m), source = 0.0d0)
        mtx(1:n, 1:n) = sys%A * sample_time
        mtx(1:n, n + 1:n + m) = sys%B * sample_time
        phi = matrix_exponential(mtx)
        ad = phi(1:n, 1:n)
        bd = phi(1:n, n + 1:n + m)
        rst = discrete_state_space(ad, bd, sys%C, sys%D, sample_time)
    case (DISCRETIZE_TUSTIN)
        ident = identity(n)
        f = mtx_inverse(ident - (0.5d0 * sample_time) * sys%A)
        ad = matmul(f, ident + (0.5d0 * sample_time) * sys%A)
        fb = matmul(f, sys%B)
        bd = sample_time * matmul(f, fb)
        dd = sys%D + (0.5d0 * sample_time) * matmul(sys%C, fb)
        rst = discrete_state_space(ad, bd, sys%C, dd, sample_time)
    case default
        error stop DYN_INVALID_INPUT_ERROR
    end select
end function

! ------------------------------------------------------------------------------
pure function d2c(sys, method) result(rst)
    !! Converts a discrete-time state-space model into a continuous-time
    !! state-space model using the inverse of the methods provided by
    !! [[c2d]].
    !! For the zero-order hold method, the continuous model is obtained from
    !! the matrix logarithm of
    !! $$ \begin{bmatrix}A_d&B_d\\0&I\end{bmatrix}. $$
    !! This requires the matrix \(A_d\) to have no eigenvalues on the
    !! negative real axis.  For the Tustin method, with \(G=(I+A_d)^{-1}\),
    !! $$ A=\tfrac{2}{T_s}G(A_d-I),\quad B=\tfrac{4}{T_s}G^2B_d,\quad
    !!    D=D_d-\tfrac{T_s}{4}C(I+A_d)B. $$
    type(discrete_state_space), intent(in) :: sys
        !! The discrete-time [[discrete_state_space]] model.
    integer, intent(in), optional :: method
        !! The discretization method, either [[DISCRETIZE_ZOH]] (the default)
        !! or [[DISCRETIZE_TUSTIN]].
    type(state_space) :: rst
        !! The resulting continuous-time [[state_space]] model.

    ! Local Variables
    integer :: m, method_, n
    real(real64) :: ts
    real(real64), allocatable, dimension(:,:) :: a, b, dc, g, ident, lg, mtx

    ! Process
    ts = sys%sample_time
    method_ = DISCRETIZE_ZOH
    if (present(method)) method_ = method
    n = size(sys%A, 1)
    m = size(sys%B, 2)
    select case (method_)
    case (DISCRETIZE_ZOH)
        allocate(mtx(n + m, n + m), source = 0.0d0)
        mtx(1:n, 1:n) = sys%A
        mtx(1:n, n + 1:n + m) = sys%B
        mtx(n + 1:n + m, n + 1:n + m) = identity(m)
        lg = matrix_log(mtx)
        a = lg(1:n, 1:n) / ts
        b = lg(1:n, n + 1:n + m) / ts
        rst = state_space(a, b, sys%C, sys%D)
    case (DISCRETIZE_TUSTIN)
        ident = identity(n)
        g = mtx_inverse(ident + sys%A)
        a = (2.0d0 / ts) * matmul(g, sys%A - ident)
        b = (4.0d0 / ts) * matmul(matmul(g, g), sys%B)
        dc = sys%D - (0.25d0 * ts) * matmul(sys%C, matmul(ident + sys%A, b))
        rst = state_space(a, b, sys%C, dc)
    case default
        error stop DYN_INVALID_INPUT_ERROR
    end select
end function

! ------------------------------------------------------------------------------
pure function z_roots(c) result(rst)
    !! Computes the roots, in the variable \(z\), of a polynomial whose
    !! coefficients are given in ascending powers of \(z^{-1}\).
    real(real64), intent(in), dimension(:) :: c
        !! The coefficients in ascending powers of \(z^{-1}\).
    complex(real64), allocatable, dimension(:) :: rst
        !! The roots of the polynomial.

    ! Local Variables
    integer :: k, n
    real(real64), allocatable, dimension(:) :: zc
    type(polynomial) :: p

    ! Process
    n = size(c)
    zc = c(n:1:-1)
    k = n
    do while (k > 1)
        if (zc(k) /= 0.0d0) exit
        k = k - 1
    end do
    if (k <= 1) then
        allocate(rst(0))
        return
    end if
    call p%initialize(zc(1:k))
    rst = p%roots()
end function

! ------------------------------------------------------------------------------
pure function matrix_sqrt(a) result(rst)
    !! Computes the principal square root of a matrix using the Denman-Beavers
    !! iteration.
    real(real64), intent(in), dimension(:,:) :: a
        !! The square matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The square root of \(A\).

    ! Local Variables
    integer :: i, n
    real(real64), allocatable, dimension(:,:) :: y, yi, ynew, z, zi

    ! Process
    n = size(a, 1)
    y = a
    z = identity(n)
    do i = 1, MAX_SQRT_ITERATIONS
        yi = mtx_inverse(y)
        zi = mtx_inverse(z)
        ynew = 0.5d0 * (y + zi)
        z = 0.5d0 * (z + yi)
        if (norm1(ynew - y) <= 1.0d-15 * norm1(ynew)) then
            y = ynew
            exit
        end if
        y = ynew
        if (i == MAX_SQRT_ITERATIONS) error stop DYN_INVALID_INPUT_ERROR
    end do
    rst = y
end function

! ------------------------------------------------------------------------------
pure function matrix_log(a) result(rst)
    !! Computes the principal matrix logarithm by taking repeated square roots
    !! until \(A^{1/2^k}\) is near the identity, then evaluating the
    !! series \( \log(I+X)=\sum_j (-1)^{j+1}X^j/j \) and scaling by
    !! \(2^k\).
    real(real64), intent(in), dimension(:,:) :: a
        !! The square matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The matrix logarithm of \(A\).

    ! Local Variables
    integer :: j, k, n
    real(real64) :: sgn
    real(real64), allocatable, dimension(:,:) :: r, term, x

    ! Process
    n = size(a, 1)
    r = a
    k = 0
    do
        x = r - identity(n)
        if (norm1(x) <= 0.25d0) exit
        if (k >= MAX_SQRT_ITERATIONS) error stop DYN_INVALID_INPUT_ERROR
        r = matrix_sqrt(r)
        k = k + 1
    end do
    allocate(rst(n, n), source = 0.0d0)
    term = identity(n)
    do j = 1, MAX_SERIES_TERMS
        term = matmul(term, x)
        sgn = merge(1.0d0, -1.0d0, mod(j, 2) == 1)
        rst = rst + (sgn / real(j, real64)) * term
        if (norm1(term) / real(j, real64) <= 1.0d-17) exit
    end do
    rst = (2.0d0 ** k) * rst
end function

! ------------------------------------------------------------------------------
pure function norm1(a) result(rst)
    !! Computes the matrix 1-norm (maximum absolute column sum).
    real(real64), intent(in), dimension(:,:) :: a
        !! The matrix.
    real(real64) :: rst
        !! The 1-norm of the matrix.

    ! Process
    rst = maxval(sum(abs(a), dim = 1))
end function

end module
