module dynamics_transfer_function
    use iso_fortran_env
    use nonlin_polynomials
    use ferror
    implicit none
    private
    public :: polynomial
    public :: transfer_function
    public :: operator(*)

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
    end type

    interface operator(*)
        module procedure :: tf_tf_mult
        module procedure :: poly_tf_mult
        module procedure :: tf_poly_mult
    end interface

contains
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