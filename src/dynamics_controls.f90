module dynamics_controls
    use iso_fortran_env
    use nonlin_polynomials
    use ferror
    use ieee_arithmetic
    use dynamics_error_handling
    use diffeq
    use linalg, only : eigen
    use lapack, only : dgetrf, dgetri, dgetrs, zgels
    implicit none
    private
    public :: polynomial
    public :: state_space
    public :: transfer_function
    public :: operator(*)
    public :: lti_solve
    public :: ss_excitation
    public :: ode_integrator

! ------------------------------------------------------------------------------
    type state_space
        !! Defines a state-space representation of a dynamic system.  This
        !! implementation takes the form:
        !!
        !! $$ \dot{x}(t) = A x(t) + B u(t) $$
        !! $$ y(t) = C x(t) + D u(t) $$
        !!
        !! Where:
        !! 
        !! - \( t \) denotes time.
        !!
        !! - \( x(t) \) is the state vector.
        !!
        !! - \( u(t) \) is the input vector.
        !!
        !! - \( y(t) \) is the output vector.
        real(real64), allocatable, dimension(:,:) :: A
            !! The N-by-N dynamics matrix, where N is the number of state
            !! variables.
        real(real64), allocatable, dimension(:,:) :: B
            !! The N-by-M input matrix, where M is the number of inputs.
        real(real64), allocatable, dimension(:,:) :: C
            !! The P-by-N output matrix, where P is the number of outputs.
        real(real64), allocatable, dimension(:,:) :: D
            !! The P-by-M feedthrough matrix.
    contains
        procedure, public :: evaluate_derivatives => ss_eval_deriv
        procedure, public :: evaluate_output => ss_eval_output
        procedure, public :: poles => ss_poles
        procedure, public :: zeros => ss_zeros
        generic, public :: transfer_function => ss_transfer_fcn, &
            ss_transfer_fcn_omega, ss_transfer_fcn_array, &
            ss_transfer_fcn_omega_array
        procedure, private :: ss_transfer_fcn
        procedure, private :: ss_transfer_fcn_omega
        procedure, private :: ss_transfer_fcn_array
        procedure, private :: ss_transfer_fcn_omega_array
    end type

    interface state_space
        module procedure :: state_space_init
        module procedure :: state_space_init_scalar
        module procedure :: state_space_init_matrices
        module procedure :: state_space_init_pid
        module procedure :: state_space_init_pid_plant
    end interface

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
        procedure, private :: tf_eval_s
        procedure, private :: tf_eval_omega
        generic, public :: evaluate => tf_eval_omega, tf_eval_s
        procedure, public :: poles => tf_poles
        procedure, public :: zeros => tf_zeros
        procedure, public :: to_ccf_state_space => tf_to_ccf_statespace
        procedure, public :: to_ocf_state_space => tf_to_ocf_statespace
    end type

    interface transfer_function
        module procedure :: init_tf_array
        module procedure :: init_tf_poly
    end interface

! ------------------------------------------------------------------------------
    interface operator(*)
        module procedure :: tf_tf_mult
        module procedure :: poly_tf_mult
        module procedure :: tf_poly_mult
        module procedure :: tf_scalar_mult
        module procedure :: scalar_tf_mult
    end interface

! ------------------------------------------------------------------------------
    interface
        subroutine ss_excitation(t, u, args)
            !! A routine for computing the excitation vector for a state-space
            !! model.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: t
                !! The time value at which to compute the excitation.
            real(real64), intent(out), dimension(:) :: u
                !! The excitation vector.
            class(*), intent(inout), optional :: args
                !! An optional argument used to pass objects in and out of the
                !! routine.
        end subroutine
    end interface

! ------------------------------------------------------------------------------
! PRIVATE TYPES
! ------------------------------------------------------------------------------
    type argument_container
        procedure(ss_excitation), pointer, nopass :: excitation
        class(*), allocatable :: user_args
        logical :: has_user_args
        class(state_space), pointer :: model
    end type

contains
! ******************************************************************************
! STATE_SPACE
! ------------------------------------------------------------------------------
pure function state_space_init(m, b, k, n_out) result(rst)
    !! Initializes the state space model.
    !!
    !! The output matrix \(C\) is initialized to one, and the
    !! feedthrough matrix \(D\) is initialized to zero.
    real(real64), intent(in), dimension(:,:) :: m
        !! The N-by-N mass matrix.
    real(real64), intent(in), dimension(size(m, 1), size(m, 2)) :: b
        !! The N-by-N damping matrix.
    real(real64), intent(in), dimension(size(m, 1), size(m, 2)) :: k
        !! The N-by-N stiffness matrix.
    integer(int32), intent(in), optional :: n_out
        !! The number of outputs.  The default is 1.
    type(state_space) :: rst
        !! The [[state_space]] model.

    ! Local Variables
    integer(int32) :: ii, jj, n, p, q, lwork, info
    integer(int32), allocatable, dimension(:) :: pvt
    real(real64) :: temp(1)
    real(real64), allocatable, dimension(:) :: work

    ! Initialization
    p = size(m, 1)
    n = 2 * p
    q = 1
    if (present(n_out)) q = n_out
    if (q <= 1) q = 1
    allocate( &
        rst%A(n, n), &
        rst%B(n, p), &
        rst%D(q, p), &
        source = 0.0d0 &
    )
    allocate(rst%C(q, n), source = 1.0d0)
    allocate(pvt(p))

    ! Workspace
    call dgetri(p, rst%b(p+1:n,:), p, pvt, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))

    ! Build p-by-p Identity sub-matrix matrix
    jj = p + 1
    do ii = 1, p
        rst%A(ii,jj) = 1.0d0
        jj = jj + 1
    end do
    
    ! Fill in the matrices
    rst%B(p+1:n,1:p) = m
    rst%A(p+1:n,1:p) = -k
    rst%A(p+1:n,p+1:n) = -b
    call dgetrf(p, p, rst%B(p+1:n,1:p), p, pvt, info)
    call dgetrs('N', p, p, rst%B(p+1:n,1:p), p, pvt, rst%A(p+1:n,1:p), p, &
        info)
    call dgetrs('N', p, p, rst%B(p+1:n,1:p), p, pvt, rst%A(p+1:n,p+1:n), &
        p, info)
    call dgetri(p, rst%B(p+1:n,1:p), p, pvt, work, lwork, info)
end function

! ------------------------------------------------------------------------------
pure function state_space_init_scalar(m, b, k) result(rst)
    !! Initializes the state space model.
    !!
    !! The output matrix \(C\) is initialized to one, and the
    !! feedthrough matrix \(D\) is initialized to zero.
    real(real64), intent(in) :: m
        !! The mass.
    real(real64), intent(in) :: b
        !! The damping.
    real(real64), intent(in) :: k
        !! The stiffness.
    type(state_space) :: rst
        !! The [[state_space]] model.

    ! Process
    allocate( &
        rst%A(2, 2), &
        rst%B(2, 1), &
        rst%D(1, 1), &
        source = 0.0d0 &
    )
    allocate(rst%C(1, 2), source = 1.0d0)
    rst%A(2,1) = -k / m
    rst%A(1,2) = 1.0d0
    rst%A(2,2) = -b / m
    rst%B(2,1) = 1.0d0 / m
end function

! ------------------------------------------------------------------------------
pure function state_space_init_matrices(a, b, c, d) result(rst)
    !! Initializes the state space model.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N dynamics matrix.
    real(real64), intent(in), dimension(:,:) :: b
        !! The N-by-M input matrix.
    real(real64), intent(in), dimension(:,:) :: c
        !! The P-by-N output matrix.
    real(real64), intent(in), dimension(:,:) :: d
        !! The P-by-M feedthrough matrix.
    type(state_space) :: rst
        !! The resulting [[state_space]] object.

    allocate(rst%A, source = a)
    allocate(rst%B, source = b)
    allocate(rst%C, source = c)
    allocate(rst%D, source = d)
end function

! ------------------------------------------------------------------------------
pure function state_space_init_pid(kp, ki, kd, tau, a, b, c, d) result(rst)
    !! Initializes a state-space model that employs a closed-loop PID
    !! controller.
    !!
    !! The PID model is augmented into the plant model as follows.
    !!
    !! $$ e = r - y $$
    !! $$ \dot{x_1} = e $$
    !! $$ \dot{x_3} = -\frac{x_3}{\tau} + \frac{e}{\tau} $$
    !! $$ u = K_{i} x_{i} + K_{p} e + \frac{K_{d}}{\tau} \left(e - 
    !! x_{3} \right) $$
    !! $$ \alpha = K_{p} + \frac{K_{d}}{\tau} $$
    !! $$ \beta = \frac{1}{1 + \alpha D} $$
    !! $$ x = \left[ \begin{matrix} x_{p} \\ x_{1} \\ x_{3} \end{matrix} 
    !! \right] $$
    !! $$ \dot{x} = A_{cl} x + B_{cl} r $$
    !! $$ y = C_{cl} x + D_{cl} r $$
    !!
    !! Where the augmented matrices are as follows.
    !!
    !! $$ A_{cl} = \left[ \begin{matrix} A - \alpha \beta B C & 
    !! \beta K_{i} B & -\frac{\beta K_{d}}{\tau} B \\ 
    !! -C + \alpha \beta D C & -\beta K_{i} D & \frac{\beta K_{d}}{\tau} D 
    !! \\ \frac{1}{\tau}\left(-C + \alpha \beta D C \right) & 
    !! -\frac{\beta K_{i}}{\tau} D & \frac{1}{\tau} \left( 
    !! 1 + \frac{\beta K_{d}}{\tau^{2}} D \right) \end{matrix} \right] $$
    !! $$ B_{cl} = \left[ \begin{matrix} \alpha \beta B \\ 1 - 
    !! \alpha \beta D \\ \frac{1}{\tau} \left( 1 - \alpha \beta D \right) 
    !! \end{matrix} \right] $$
    !! $$ C_{cl} = \left[ \begin{matrix} C - \alpha \beta D C & 
    !! \beta K_{i} D & -\frac{\beta K_{d}}{\tau} D \end{matrix} \right] $$
    !! $$ D_{cl} = \alpha \beta D $$
    real(real64), intent(in) :: kp
        !! The proportional gain term.
    real(real64), intent(in) :: ki
        !! The integral gain term.
    real(real64), intent(in) :: kd
        !! The derivative gain term.
    real(real64), intent(in) :: tau
        !! The time constant of the first order derivative filter
        !! \( K_{d} \frac{s}{\tau s + 1} \).
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N dynamics matrix for the plant.
    real(real64), intent(in), dimension(size(a, 1), 1) :: b
        !! The N-by-1 input matrix for the plant.
    real(real64), intent(in), dimension(1, size(a, 1)) :: c
        !! The 1-by-N output matrix for the plant.
    real(real64), intent(in), dimension(1, 1) :: d
        !! The 1-by-1 feedthrough matrix for the plant.
    type(state_space) :: rst
        !! The resulting [[state_space]] object.

    ! Local Variables
    integer(int32) :: n, n1, n2
    real(real64) :: alpha, beta, ab
    real(real64), allocatable, dimension(:,:) :: dc

    ! Initialization
    n = size(a, 1)
    n1 = n + 1
    n2 = n + 2
    allocate( &
        rst%A(n2, n2), &
        rst%B(n2, 1), &
        rst%C(1, n2), &
        rst%D(1, 1) &
    )
    dc = matmul(d, c)

    ! Process
    alpha = kp + kd / tau
    if (d(1,1) == 0.0d0) then
        beta = 1.0d0
    else
        beta = 1.0d0 / (1.0d0 + alpha * d(1,1))
    end if
    ab = alpha * beta
    rst%A(1:n,1:n) = a - ab * matmul(b, c)
    rst%A(n1,1:n) = -c(1,1:n) + ab * dc(1,1:n)
    rst%A(4,1:n) = (1.0d0 / tau) * rst%A(3,1:2)

    rst%A(1:n,n1) = beta * ki * b(1:n,1)
    rst%A(n1,n1) = -beta * ki * d(1,1)
    rst%A(n2,n1) = rst%A(3,3) / tau

    rst%A(1:n,n2) = -(beta * kd / tau) * b(1:n,1)
    rst%A(n1,n2) = beta * kd * d(1,1) / tau
    rst%A(n2,n2) = -1.0d0 / tau + beta * kd * d(1,1) / (tau**2)

    rst%B(1:n,1) = ab * b(1:n,1)
    rst%B(n1,1) = 1.0d0 - ab * d(1,1)
    rst%B(n2,1) = (1.0d0 / tau) * rst%B(3,1)

    rst%C(1,1:n) = c(1,1:n) - ab * dc(1,1:n)
    rst%C(1,n1) = beta * ki * d(1,1)
    rst%C(1,n2) = -(beta * kd / tau) * d(1,1)

    rst%D = ab * d
end function

! ------------------------------------------------------------------------------
pure function state_space_init_pid_plant(kp, ki, kd, tau, plant) result(rst)
    !! Initializes a state-space model that employs a closed-loop PID
    !! controller.
    !!
    !! The PID model is augmented into the plant model as follows.
    !!
    !! $$ e = r - y $$
    !! $$ \dot{x_1} = e $$
    !! $$ \dot{x_3} = -\frac{x_3}{\tau} + \frac{e}{\tau} $$
    !! $$ u = K_{i} x_{i} + K_{p} e + \frac{K_{d}}{\tau} \left(e - 
    !! x_{3} \right) $$
    !! $$ \alpha = K_{p} + \frac{K_{d}}{\tau} $$
    !! $$ \beta = \frac{1}{1 + \alpha D} $$
    !! $$ x = \left[ \begin{matrix} x_{p} \\ x_{1} \\ x_{3} \end{matrix} 
    !! \right] $$
    !! $$ \dot{x} = A_{cl} x + B_{cl} r $$
    !! $$ y = C_{cl} x + D_{cl} r $$
    !!
    !! Where the augmented matrices are as follows.
    !!
    !! $$ A_{cl} = \left[ \begin{matrix} A - \alpha \beta B C & 
    !! \beta K_{i} B & -\frac{\beta K_{d}}{\tau} B \\ 
    !! -C + \alpha \beta D C & -\beta K_{i} D & \frac{\beta K_{d}}{\tau} D 
    !! \\ \frac{1}{\tau}\left(-C + \alpha \beta D C \right) & 
    !! -\frac{\beta K_{i}}{\tau} D & \frac{1}{\tau} \left( 
    !! 1 + \frac{\beta K_{d}}{\tau^{2}} D \right) \end{matrix} \right] $$
    !! $$ B_{cl} = \left[ \begin{matrix} \alpha \beta B \\ 1 - 
    !! \alpha \beta D \\ \frac{1}{\tau} \left( 1 - \alpha \beta D \right) 
    !! \end{matrix} \right] $$
    !! $$ C_{cl} = \left[ \begin{matrix} C - \alpha \beta D C & 
    !! \beta K_{i} D & -\frac{\beta K_{d}}{\tau} D \end{matrix} \right] $$
    !! $$ D_{cl} = \alpha \beta D $$
    real(real64), intent(in) :: kp
        !! The proportional gain term.
    real(real64), intent(in) :: ki
        !! The integral gain term.
    real(real64), intent(in) :: kd
        !! The derivative gain term.
    real(real64), intent(in) :: tau
        !! The time constant of the first order derivative filter
        !! \( K_{d} \frac{s}{\tau s + 1} \).
    class(state_space), intent(in) :: plant
        !! The plant model.
    type(state_space) :: rst
        !! The resulting [[state_space]] object.

    rst = state_space_init_pid(kp, ki, kd, tau, plant%A, plant%B, &
        plant%C, plant%D)
end function

! ------------------------------------------------------------------------------
pure function ss_eval_deriv(this, u, x) result(rst)
    !! Evaluates the state time derivative \( \dot{x}(t) = A x(t) + 
    !! B u(t) \).
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    real(real64), intent(in), dimension(:) :: u
        !! The M-element input array.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element state array.
    real(real64), allocatable, dimension(:) :: rst
        !! The N-element state time derivative vector.

    rst = matmul(this%A, x) + matmul(this%B, u)
end function

! ------------------------------------------------------------------------------
pure function ss_eval_output(this, u, x) result(rst)
    !! Evaluates the output vector \( y(t) = C x(t) + D u(t) \).
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    real(real64), intent(in), dimension(:) :: u
        !! The M-element input array.
    real(real64), intent(in), dimension(:) :: x
        !! The N-element state array.
    real(real64), allocatable, dimension(:) :: rst
        !! The P-element output array.

    rst = matmul(this%C, x) + matmul(this%D, u)
end function

! ------------------------------------------------------------------------------
function ss_poles(this, err) result(rst)
    !! Computes the poles of the state space model.
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The poles of the model.

    ! Local Variables
    integer(int32) :: n
    real(real64), allocatable, dimension(:,:) :: ac

    ! Process
    n = size(this%A, 1)
    allocate(rst(n))
    if (n == 0) return
    allocate(ac(n, n), source = this%A)
    call eigen(ac, rst, err = err)
end function

! ------------------------------------------------------------------------------
function ss_zeros(this, err) result(rst)
    !! Computes the zeros of the state space model.
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    complex(real64), allocatable, dimension(:) :: rst
        !! The zeros of the model.

    ! Local Variables
    integer(int32) :: i, j, n, m, p, nz
    real(real64), allocatable, dimension(:,:) :: Az, Bz
    complex(real64), allocatable, dimension(:) :: buffer, trimmed

    ! Initialization
    n = size(this%A, 1)
    m = size(this%B, 2)
    p = size(this%C, 1)
    nz = n + max(m, p)
    allocate(Az(nz, nz), Bz(nz, nz), source = 0.0d0)
    allocate(buffer(nz), trimmed(nz))

    ! Build the matrices
    Az(1:n,1:n) = this%A
    Az(1:n,n+1:n+m) = this%B
    do i = 1, nz
        Bz(i,i) = 1.0d0
    end do
    Bz(n+1:n+p,1:n) = this%C
    Bz(n+1:n+p,n+1:n+m) = this%D

    ! Compute the eigenvalues
    call eigen(Az, Bz, buffer)

    ! Only keep the physically realizable zeros - exclude NaN's and (0,0)
    j = 0
    do i = 1, nz
        if (ieee_is_nan(real(buffer(i)))) cycle
        if (real(buffer(i)) == 0.0d0 .and. aimag(buffer(i)) == 0.0d0) cycle
        j = j + 1
        trimmed(j) = buffer(i)
    end do
    rst = trimmed(1:j)
end function

! ------------------------------------------------------------------------------
pure function ss_transfer_fcn(this, s) result(rst)
    !! Evaluates the transfer functions for the model at the parameter \(s\).
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    complex(real64), intent(in) :: s
        !! The frequency at which to evaluate the transfer functions.
    complex(real64), allocatable, dimension(:,:) :: rst
        !! The resulting transfer functions.

    ! Local Variables
    integer(int32) :: j, ninput, noutput, n, info, lwork
    complex(real64) :: temp(1)
    complex(real64), allocatable, dimension(:) :: work, Bj
    complex(real64), allocatable, dimension(:,:) :: A, Ac

    ! Initialization
    n = size(this%A, 1)
    ninput = size(this%B, 2)
    noutput = size(this%C, 1)
    allocate(rst(noutput, ninput))
    allocate(Ac(n, n), Bj(n))

    ! Determine DGELS workspace requirements
    call zgels('N', n, n, 1, Ac, n, Bj, n, temp, -1, info)
    lwork = int(temp(1), int32)
    allocate(work(lwork))

    ! Build s * I - A
    do j = 1, n
        Ac(:,j) = -this%A(:,j)
        Ac(j,j) = Ac(j,j) + s
    end do
    allocate(A(n, n), source = Ac)

    ! For each column j, solve Ac * x(j) = B(j) and then compute
    ! the output: C x(j) + D(j)
    do j = 1, ninput
        if (j /= 1) A = Ac  ! We need a copy as A will be overwritten
        Bj = this%B(:,j)
        call zgels('N', n, n, 1, A, n, Bj, n, work, lwork, info)
        rst(:,j) = matmul(this%C, Bj) + this%D(:,j)
    end do
end function
    
! ------------------------------------------------------------------------------
pure function ss_transfer_fcn_omega(this, omega) result(rst)
    !! Evaluates the transfer functions for the model at frequency \(\omega\).
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    real(real64), intent(in) :: omega
        !! The frequency at which to evaluate the transfer functions.
    complex(real64), allocatable, dimension(:,:) :: rst
        !! The resulting transfer functions.

    ! Local Variables
    complex(real64), parameter :: j = (0.0d0, 1.0d0)
    complex(real64) :: s

    ! Process
    s = j * omega
    rst = ss_transfer_fcn(this, s)
end function

! ------------------------------------------------------------------------------
pure function ss_transfer_fcn_array(this, s) result(rst)
    !! Evaluates the transfer functions for the model at the frequencies given
    !! in the array \(s\).
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    complex(real64), intent(in), dimension(:) :: s
        !! The frequencies at which to evaluate the transfer functions.
    complex(real64), allocatable, dimension(:,:,:) :: rst
        !! The resulting transfer functions, with each page of the array 
        !! containing the transfer functions for a specific frequency.

    ! Local Variables
    integer(int32) :: i, ninput, noutput, n

    ! Initialization
    ninput = size(this%B, 2)
    noutput = size(this%C, 1)
    n = size(s)
    allocate(rst(noutput, ninput, n))

    ! Process
    do concurrent (i = 1:n)
        rst(:,:,i) = ss_transfer_fcn(this, s(i))
    end do
end function

! ------------------------------------------------------------------------------
pure function ss_transfer_fcn_omega_array(this, omega) result(rst)
    !! Evaluates the transfer functions for the model at the frequencies given
    !! in the array \(omega\).
    class(state_space), intent(in) :: this
        !! The [[state_space]] object.
    real(real64), intent(in), dimension(:) :: omega
        !! The frequencies at which to evaluate the transfer functions.
    complex(real64), allocatable, dimension(:,:,:) :: rst
        !! The resulting transfer functions, with each page of the array 
        !! containing the transfer functions for a specific frequency.

    ! Local Variables
    integer(int32) :: i, ninput, noutput, n

    ! Initialization
    ninput = size(this%B, 2)
    noutput = size(this%C, 1)
    n = size(omega)
    allocate(rst(noutput, ninput, n))

    ! Process
    do concurrent (i = 1:n)
        rst(:,:,i) = ss_transfer_fcn_omega(this, omega(i))
    end do
end function

! ******************************************************************************
! TRANSFER_FUNCTION
! ------------------------------------------------------------------------------
function init_tf_poly(y, x) result(rst)
    !! Initializes a new transfer function.
    class(polynomial), intent(in) :: y
        !! The numerator polynomial \(Y(s)\) in 
        !! \( H(s) = \frac{Y(s)}{X(s)} \).
    class(polynomial), intent(in) :: x
        !! The denominator polynomial \(X(s)\) in 
        !! \( H(s) = \frac{Y(s)}{X(s)} \).
    type(transfer_function) :: rst
        !! The resulting [[transfer_function]].

    ! Process
    call rst%Y%initialize(y%get_all())
    call rst%X%initialize(x%get_all())
end function

! ------------------------------------------------------------------------------
function init_tf_array(y, x) result(rst)
    !! Initializes a new transfer function.
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
    type(transfer_function) :: rst
        !! The resulting [[transfer_function]].

    ! Process
    call rst%Y%initialize(y)
    call rst%X%initialize(x)
end function

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
    !! \(\omega\).
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
function tf_scalar_mult(x, y) result(rst)
    !! Multiplies a transfer function by a scalar value.
    class(transfer_function), intent(in) :: x
        !! The left-hand-side argument.
    real(real64), intent(in) :: y
        !! The right-hand-side argument.
    type(transfer_function) :: rst
        !! The resulting transfer function.

    ! Process
    rst%Y = y * x%Y
    rst%X = x%X
end function

! ------------------------------------------------------------------------------
function scalar_tf_mult(x, y) result(rst)
    !! Multiplies a transfer function by a scalar value.
    real(real64), intent(in) :: x
        !! The left-hand-side argument.
    class(transfer_function), intent(in) :: y
        !! The right-hand-side argument.
    type(transfer_function) :: rst
        !! The resulting transfer function.

    ! Process
    rst%Y = x * y%Y
    rst%X = y%X
end function

! ******************************************************************************
! LTI SOLVERS
! ------------------------------------------------------------------------------
function lti_solve(mdl, u, t, ic, solver, args, err) result(rst)
    !! Solves the LTI system given by the specified state space model.
    class(state_space), intent(in), target :: mdl
        !! The state_space model to solve.
    procedure(ss_excitation), pointer, intent(in) :: u
        !! The routine used to compute the excitation vector.
    real(real64), intent(in), dimension(:) :: t
        !! The time points at which to compute the solution.  The array must
        !! have at least 2 values; however, more may be specified.  If only
        !! 2 values are specified, the integrator will compute the solution at
        !! those points, but it will also return any intermediate integration
        !! steps that may be required.  However, if more than 2 points are
        !! given, the integrator will return the solution values only at the
        !! specified time points.
    real(real64), intent(in), dimension(:) :: ic
        !! The initial condition vector.  This array must be the same size as
        !! the number of state variables.
    class(ode_integrator), intent(in), optional, target :: solver
        !! The ODE solver to utilize.  If not specified, the default solver
        !! is a 4th/5th order Runge-Kutta integrator.
    class(*), intent(inout), optional :: args
        !! An optional container for arguments to pass to the excitation
        !! routine.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The solution.  The time points at which the solution was evaluated
        !! are stored in the first column and the output(s) are stored in the
        !! remaining column(s).

    ! Local Variables
    integer(int32) :: i, n, npts, flag, nInputs, nOutputs
    real(real64), allocatable, dimension(:) :: uv
    real(real64), allocatable, dimension(:,:) :: sol
    type(argument_container) :: container
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    class(ode_integrator), pointer :: integrator
    type(runge_kutta_45), target :: defaultIntegrator
    type(ode_container) :: odeMdl
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    container%excitation => u
    container%model => mdl
    container%has_user_args = present(args)
    if (present(args)) allocate(container%user_args, source = args)
    n = size(mdl%A, 1)
    nInputs = size(mdl%B, 2)
    nOutputs = size(mdl%C, 1)

    ! Input Checking
    if (size(ic) /= n) then
        call report_array_size_error("lti_solve_statespace", "ic", n, &
            size(ic), errmgr)
        return
    end if
    
    ! Set up the integrator
    if (present(solver)) then
        integrator => solver
    else
        integrator => defaultIntegrator
    end if
    odeMdl%fcn => ode_solver_routine

    ! Call the solver
    call integrator%solve(odeMdl, t, ic, args = container, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Get the output from the solver
    sol = integrator%get_solution()
    npts = size(sol, 1)
    allocate(rst(npts, nOutputs + 1), uv(nInputs), stat = flag, source = 0.0d0)
    if (flag /= 0) then
        call report_memory_error("lti_solve_statespace", flag, errmgr)
        return
    end if

    ! Compute: C * x + D * u
    do i = 1, npts
        call u(sol(i,1), uv, args)
        rst(i,1) = sol(i,1)
        rst(i,2:) = mdl%evaluate_output(uv, sol(i,2:))
    end do
end function

! ----------
subroutine ode_solver_routine(t, x, dxdt, args)
    ! The routine called by the LTI solver
    real(real64), intent(in) :: t, x(:)
    real(real64), intent(out) :: dxdt(:)
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: n, p
    real(real64) :: nan
    real(real64), allocatable, dimension(:) :: u

    ! Initialization
    nan = ieee_value(nan, IEEE_QUIET_NAN)
    
    ! Process
    select type (args)
    class is (argument_container)
        ! Evaluate the excitation vector at t
        n = size(args%model%B, 2)
        allocate(u(n), source = 0.0d0)
        if (args%has_user_args) then
            call args%excitation(t, u, args = args%user_args)
        else
            call args%excitation(t, u)
        end if

        ! Evaluate the model
        dxdt = args%model%evaluate_derivatives(u, x)
    class default
        dxdt = nan
    end select
end subroutine

! ------------------------------------------------------------------------------
end module