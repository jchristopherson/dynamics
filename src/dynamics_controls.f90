module dynamics_controls
    use iso_fortran_env
    use nonlin_polynomials
    use ferror
    implicit none
    private
    public :: polynomial
    public :: state_space
    public :: transfer_function
    public :: operator(*)

! ------------------------------------------------------------------------------
    type state_space
        !! Defines a state-space representation of a dynamic system.  This
        !! implementation takes the form:
        !!
        !! $$ \dot{x}(t) = A(t) x(t) + B(t) u(t) $$
        !! $$ y(t) = C(t) x(t) + D(t) u(t) $$
        !!
        !! Where:
        !! 
        !! - /( t /) denotes time.
        !!
        !! - /( x(t) /) is the state vector.
        !!
        !! - /( u(t) /) is the input vector.
        !!
        !! - /( y(t) /) is the output vector.
        real(real64), allocatable, dimension(:,:) :: A
            !! The N-by-N dynamics matrix, where N is the number of state
            !! variables.
        real(real64), allocatable, dimension(:,:) :: B
            !! The N-by-M input matrix, where M is the number of inputs.
        real(real64), allocatable, dimension(:,:) :: C
            !! The P-by-N output matrix, where P is the number of outputs.
        real(real64), allocatable, dimension(:,:) :: D
            !! The P-by-M feedthrough matrix.
    end type

! ------------------------------------------------------------------------------
    type transfer_function
        !! Defines a transfer function for a continuous system of the form
        !! \( H(s) = \frac{Y(s)}{X(s)} \).
        type(polynomial) :: Y
            !! The numerator polynomial \(Y(s)\) in 
            !! \( H(s) = \frac{Y(s)}{X(s)} \).  The polynomial coefficients
            !! are stored in acending order such that 
            !! \( y_1 + y_2 s + y_3 s^2 ... \).
        type(polynomial) :: X
            !! The denominator polynomial \(X(s)\) in 
            !! \( H(s) = \frac{Y(s)}{X(s)} \).  The polynomial coefficients
            !! are stored in acending order such that 
            !! \( x_1 + x_2 s + x_3 s^2 ... \).
    contains
        procedure, private :: tf_init_poly
        procedure, private :: tf_init_array
        generic, public :: initialize => tf_init_poly, tf_init_array
        procedure, private :: tf_eval_s
        procedure, private :: tf_eval_omega
        generic, public :: evaluate => tf_eval_omega, tf_eval_s
        procedure, public :: poles => tf_poles
        procedure, public :: zeros => tf_zeros
        procedure, public :: to_ccf_state_space => tf_to_ccf_statespace
        procedure, public :: to_ocf_state_space => tf_to_ocf_statespace
    end type

! ------------------------------------------------------------------------------
    interface operator(*)
        module procedure :: tf_tf_mult
        module procedure :: poly_tf_mult
        module procedure :: tf_poly_mult
    end interface

contains
! ******************************************************************************
! TRANSFER_FUNCTION
! ------------------------------------------------------------------------------
subroutine tf_init_poly(this, y, x)
    !! Initializes a new transfer function.
    class(transfer_function), intent(inout) :: this
        !! The transfer_function object.
    class(polynomial), intent(in) :: y
        !! The numerator polynomial \(Y(s)\) in 
        !! \( H(s) = \frac{Y(s)}{X(s)} \).
    class(polynomial), intent(in) :: x
        !! The denominator polynomial \(X(s)\) in 
        !! \( H(s) = \frac{Y(s)}{X(s)} \).

    ! Process
    call this%Y%initialize(y%get_all())
    call this%X%initialize(x%get_all())
end subroutine

! ------------------------------------------------------------------------------
subroutine tf_init_array(this, y, x)
    !! Initializes a new transfer function.
    class(transfer_function), intent(inout) :: this
        !! The transfer_function object.
    real(real64), intent(in), dimension(:) :: y
        !! The numerator polynomial \(Y(s)\) in 
        !! \( H(s) = \frac{Y(s)}{X(s)} \).  The polynomial coefficients
        !! are stored in acending order such that 
        !! \( y_1 + y_2 s + y_3 s^2 ... \).
    real(real64), intent(in), dimension(:) :: x
        !! The denominator polynomial \(X(s)\) in 
        !! \( H(s) = \frac{Y(s)}{X(s)} \).  The polynomial coefficients
        !! are stored in acending order such that 
        !! \( x_1 + x_2 s + x_3 s^2 ... \).

    ! Process
    call this%Y%initialize(y)
    call this%X%initialize(x)
end subroutine

! ------------------------------------------------------------------------------
pure elemental function tf_eval_s(this, s) result(rst)
    !! Evaluates the transfer function at the specified value of the
    !! Laplace variable \(s\).
    class(transfer_function), intent(in) :: this
        !! The transfer_function object.
    complex(real64), intent(in) :: s
        !! The Laplace variable at which to evaluate the transfer function.
    complex(real64) :: rst
        !! The value of the transfer function.

    ! Process
    rst = this%Y%evaluate(s) / this%X%evaluate(s)
end function

! ------------------------------------------------------------------------------
pure elemental function tf_eval_omega(this, omega) result(rst)
    !! Evaluates the transfer function at the specified value of the
    !! Laplace variable \(s\).
    class(transfer_function), intent(in) :: this
        !! The transfer_function object.
    real(real64), intent(in) :: omega
        !! The frequency, in rad/s, at which to evaluate the transfer 
        !! function.
    complex(real64) :: rst
        !! The value of the transfer function.

    ! Parameters
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    complex(real64) :: s

    ! Process
    s = j * omega
    rst = this%tf_eval_s(s)
end function

! ------------------------------------------------------------------------------
function tf_poles(this, err) result(rst)
    !! Computes the poles of the transfer function.
    class(transfer_function), intent(in) :: this
        !! The transfer_function object.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The poles of the transfer function.

    ! Process
    rst = this%X%roots(err = err)
end function

! ------------------------------------------------------------------------------
function tf_zeros(this, err) result(rst)
    !! Computes the zeros of the transfer function.
    class(transfer_function), intent(in) :: this
        !! The transfer function object.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The zeros of the transfer function.

    ! Process
    rst = this%Y%roots(err = err)
end function

! ------------------------------------------------------------------------------
function tf_to_ccf_statespace(this) result(rst)
    !! Converts a transfer_function type into a controllable canonical form
    !! state_space type.  See 
    !! [this](https://en.wikipedia.org/wiki/State-space_representation) article
    !! for a description of this form.
    class(transfer_function), intent(in) :: this
        !! The transfer_function to convert.
    type(state_space) :: rst
        !! The resulting state-space object.

    ! Local Variables
    integer(int32) :: i, order, n, norder
    real(real64) :: a

    ! Process
    order = this%X%order()
    n = order + 1
    allocate(rst%A(order, order), rst%B(order, 1), rst%C(1, order), &
        rst%D(1, 1), source = 0.0d0)
    a = this%X%get(n)
    do i = 1, order
        rst%A(order,i) = -this%X%get(i) / a
    end do
    do i = 1, order - 1
        rst%A(i,i+1) = 1.0d0
    end do
    rst%B(order, 1) = 1.0d0

    norder = this%Y%order()
    do i = 1, min(norder + 1, order)
        rst%C(1,i) = this%Y%get(i) / a
    end do
end function

! ------------------------------------------------------------------------------
function tf_to_ocf_statespace(this) result(rst)
    !! Converts a transfer_function type into an observable canonical form
    !! state_space type.  See 
    !! [this](https://en.wikipedia.org/wiki/State-space_representation) article
    !! for a description of this form.
    class(transfer_function), intent(in) :: this
        !! The transfer_function to convert.
    type(state_space) :: rst
        !! The resulting state-space object.

    ! Local Variables
    integer(int32) :: i, order, n, norder
    real(real64) :: a

    ! Process
    order = this%X%order()
    n = order + 1
    allocate(rst%A(order, order), rst%B(order, 1), rst%C(1, order), &
        rst%D(1, 1), source = 0.0d0)
    a = this%X%get(n)
    do i = 1, order
        rst%A(i,order) = -this%X%get(i) / a
    end do
    do i = 1, order - 1
        rst%A(i+1,i) = 1.0d0
    end do

    norder = this%Y%order()
    do i = 1, min(norder + 1, order)
        rst%B(i,1) = this%Y%get(i) / a
    end do

    rst%C(1,order) = 1.0d0
end function

! ******************************************************************************
! OPERATORS
! ------------------------------------------------------------------------------
function tf_tf_mult(x, y) result(rst)
    !! Multiplies two transfer functions.
    class(transfer_function), intent(in) :: x
        !! The left-hand-side argument.
    class(transfer_function), intent(in) :: y
        !! The right-hand-side argument.
    type(transfer_function) :: rst
        !! The resulting transfer function.

    ! Process
    rst%Y = x%Y * y%Y
    rst%X = x%X * y%X
end function

! ------------------------------------------------------------------------------
function poly_tf_mult(x, y) result(rst)
    !! Multiplies a polynomial and a transfer function to result in a new
    !! transfer function.
    class(polynomial), intent(in) :: x
        !! The left-hand-side argument.
    class(transfer_function), intent(in) :: y
        !! The right-hand-side argument.
    type(transfer_function) :: rst
        !! The resulting transfer function.
    
    ! Process
    rst%Y = x * y%Y
    rst%X = y%X
end function

! ------------------------------------------------------------------------------
function tf_poly_mult(x, y) result(rst)
    !! Multiplies a transfer function and a polynomial to result in a new
    !! transfer function.
    class(transfer_function), intent(in) :: x
        !! The left-hand-side argument.
    class(polynomial), intent(in) :: y
        !! The right-hand-side argument.
    type(transfer_function) :: rst
        !! The resulting transfer function.

    ! Process
    rst%Y = x%Y * y
    rst%X = x%X
end function

! ------------------------------------------------------------------------------
end module