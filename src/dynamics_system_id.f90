module dynamics_system_id
    use iso_fortran_env
    use ferror
    use dynamics_error_handling
    use diffeq
    use fstats
    implicit none
    private
    public :: dynamic_system_measurement
    public :: model_information
    public :: siso_model_fit_least_squares
    public :: ode_integrator
    public :: regression_statistics
    public :: iteration_controls
    public :: lm_solver_options
    public :: convergence_info
    public :: iteration_update
    public :: constraint_equations


    interface
        subroutine constraint_equations(xg, fg, xc, p, fc, args)
            !! An interface to a set of routines for defining constraint 
            !! equations to the fitting process.
            use iso_fortran_env, only : real64
            real(real64), intent(in), dimension(:) :: xg
                !! An N-element array containing the N independent variable
                !! values for the N differential equation solution points.
            real(real64), intent(in), dimension(:) :: fg
                !! An N-element array containing the N differential equation
                !! solution points.
            real(real64), intent(in), dimension(:) :: xc
                !! An M-element array containing the M independent variable
                !! values for the M constraint equations.
            real(real64), intent(in), dimension(:) :: p
                !! An array containing the model parameters.
            real(real64), intent(out), dimension(:) :: fc
                !! An M-element array where the values of the constraint 
                !! equations should be written.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass data in/out
                !! of this routine.
        end subroutine
    end interface

    type dynamic_system_measurement
        !! A container of a single measurement data set.
        real(real64), allocatable, dimension(:) :: t
            !! The time points at which the measurements were taken.
        real(real64), allocatable, dimension(:) :: output
            !! The output data.
        real(real64), allocatable, dimension(:) :: input
            !! The input data.
    end type

    type model_information
        !! A container for model information.
        real(real64), allocatable, dimension(:) :: model
            !! An array containing the model parameters.
        class(base_interpolator), pointer :: excitation
            !! An interpolation object allowing sampling of the excitation 
            !! function.
        class(*), pointer :: user_info
            !! Information the user has passed along.
    end type

    type regression_information
        !! A container of information being shared with the regression function.
        integer(int32), allocatable, dimension(:,:) :: start_stop
            !! An N-by-2 matrix containing the starting indices of each data
            !! set in the first column, and the ending indices of each data
            !! set in the second column.
        class(ode_integrator), pointer :: integrator
            !! A pointer to the ODE integrator.
        integer(int32) :: solution_index
            !! The column index of the solution component to extract from the
            !! ODE solver output.  Do not account for the time vector (e.g.
            !! enter a value of 1 to extract the first solution component.)
        real(real64), allocatable, dimension(:,:) :: initial_conditions
            !! The initial conditions array to pass to the ODE solver.
        real(real64), allocatable, dimension(:) :: excitation_data
            !! An array of excitation data to pass to the ODE solver.
        procedure(ode), pointer, nopass :: ode_routine
            !! The routine, defined by the calling code, containing the ODE's
            !! to solve.
        logical :: uses_constraints
            !! True if constraint equations are to be utilized; else, false.
        procedure(constraint_equations), pointer, nopass :: constraints
            !! A pointer to the constraint equations routine.
        class(*), pointer :: user_info
            !! Information the user has passed along.
    end type

    interface siso_model_fit_least_squares
        module procedure :: siso_model_fit_least_squares_1
        module procedure :: siso_model_fit_least_squares_2
    end interface

contains
! ------------------------------------------------------------------------------
subroutine siso_model_fit_least_squares_1(fcn, x, ic, p, integrator, ind, &
    maxp, minp, stats, alpha, controls, settings, info, status, cov, xc, yc, &
    constraints, weights, args, &
    err)
    !! Attempts to fit a model of a single-intput, single-output (SISO) dynamic 
    !! system by means of an iterative least-squares solver.  The algorithm
    !! computes the solution to the differential equations numerically, and
    !! compares the output to the known solution via a Levenberg-Marquardt
    !! least-squares solver.
    procedure(ode), pointer, intent(in) :: fcn
        !! The routine containing the ODE's being fit.  To communicate 
        !! model parameters and other relevant information, an instance of the
        !! [[model_information]] type is passed to the optional argument of this
        !! routine.  Use the "select type" construct to access this information.
    class(dynamic_system_measurement), intent(in), dimension(:) :: x
        !! An M-element array of arrays with each array containing the measured
        !! input and output of the system being identified.
    real(real64), intent(in), dimension(:) :: ic
        !! The initial condition vector for the equations in fcn.
    real(real64), intent(inout), dimension(:) :: p
        !! An N-element array containing an initial guess at the parameters.  
        !! On output, the computed model parameters.
    class(ode_integrator), intent(inout), optional, target :: integrator
        !! The integrator to use when solving the system equations.  If not
        !! supplied, the default integrator will be used.  The default 
        !! integrator is a Runge-Kutta integrator (Dormand-Prince).
    integer(int32), intent(in), optional :: ind
        !! The index of the ODE in fcn providing the output to fit.  If
        !! no value is supplied, a value of 1 will be utilized.
    real(real64), intent(in), optional, dimension(:) :: maxp
        !! An optional N-element array that can be used as upper limits on the 
        !! parameter values. If no upper limit is requested for a particular 
        !! parameter, utilize a very large value. The internal default is to 
        !! utilize huge() as a value.
    real(real64), intent(in), optional, dimension(:) :: minp
        !! An optional N-element array that can be used as lower limits on the 
        !! parameter values. If no lower limit is requested for a particalar 
        !! parameter, utilize a very large magnitude, but negative, value. The 
        !! internal default is to utilize -huge() as a value.
    type(regression_statistics), intent(out), optional, dimension(:) :: stats
        !! An optional N-element array that, if supplied, will be used to 
        !! return statistics about the fit for each parameter.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals. The default value is 0.05 such that a 95% confidence 
        !! interval is calculated.
    type(iteration_controls), intent(in), optional :: controls
        !! An optional input providing custom iteration controls.
    type(lm_solver_options), intent(in), optional :: settings
        !! An optional input providing custom settings for the solver.
    type(convergence_info), intent(out), optional :: info
        !! An optional output that can be used to gain information about the 
        !! iterative solution and the nature of the convergence.
    procedure(iteration_update), pointer, intent(in), optional :: status
        !! An optional pointer to a routine that can be used to extract 
        !! iteration information.
    real(real64), intent(out), optional, dimension(:,:) :: cov
        !! An optional N-by-N matrix that, if supplied, will be used to return 
        !! the covariance matrix.
    real(real64), intent(in), optional, dimension(:) :: xc
        !! An optional NC-element array containing the values of the independent 
        !! variable at which the constraint equations are defined.
    real(real64), intent(in), optional, dimension(:) :: yc
        !! An optional NC-element array containing the constraint function 
        !! values at xc.
    procedure(constraint_equations), pointer, optional :: constraints
        !! An optional input, that must be utilized with the xc and yc inputs,
        !! but allows for the implementation of additional constraints on the
        !! solution outside of the differential equations being fitted.  An
        !! example usage would be an additional set of quasi-static tests that
        !! could help identify a stiffness term, for instance.  Other uses of
        !! course can be imagined.
    real(real64), intent(in), optional, dimension(:) :: weights
        !! An optional array containing weighting factors for every equation.
    class(*), intent(inout), optional, target :: args
        !! User-defined information to pass along to fcn.  These arguments,
        !! if supplied, will be passed through to fcn by means of the
        !! [[model_information]] type.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, i1, i2, n, ni, npts, nc, flag
    real(real64), allocatable, dimension(:) :: t, f, ymod, resid
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    type(runge_kutta_45), target :: default_ode_solver
    type(regression_information) :: addinfo
    procedure(regression_function), pointer :: fcnptr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(xc) .and. present(yc) .and. present(constraints)) then
        nc = size(xc)
        if (size(yc) /= nc) then
            call report_array_size_error("siso_model_fit_least_squares_1", &
                "yc", nc, size(yc), errmgr)
            return
        end if
    else
        nc = 0
    end if
    n = size(x)
    fcnptr => nlsq_fun

    allocate(addinfo%start_stop(n, 2), stat = flag)
    if (flag /= 0) go to 100
    i1 = 1
    i2 = 0
    npts = 0
    do i = 1, n
        ni = size(x(i)%t)
        npts = npts + ni
        i2 = i2 + ni
        addinfo%start_stop(i, 1) = i1
        addinfo%start_stop(i, 2) = i2
        i1 = i2 + 1
    end do
    npts = npts + nc
    allocate(t(npts), f(npts), addinfo%excitation_data(npts), ymod(npts), &
        resid(npts), stat = flag)
    if (flag /= 0) go to 100

    do i = 1, n
        i1 = addinfo%start_stop(i, 1)
        i2 = addinfo%start_stop(i, 2)
        t(i1:i2) = x(i)%t
        addinfo%excitation_data(i1:i2) = x(i)%input
        f(i1:i2) = x(i)%output
    end do
    if (nc > 0) then
        i1 = i2 + 1
        t(i1:npts) = xc
        f(i1:npts) = yc
    end if

    if (present(integrator)) then
        addinfo%integrator => integrator
    else
        addinfo%integrator => default_ode_solver
    end if

    allocate(addinfo%initial_conditions(1, size(ic)))
    addinfo%initial_conditions(1,:) = ic

    if (present(ind)) then
        addinfo%solution_index = ind
    else
        addinfo%solution_index = 1
    end if

    addinfo%ode_routine => fcn

    if (nc > 0) then
        addinfo%uses_constraints = .true.
        addinfo%constraints => constraints
    else
        addinfo%uses_constraints = .false.
    end if

    if (present(args)) addinfo%user_info => args

    ! Solve the system
    call nonlinear_least_squares(fcnptr, t, f, p, ymod, resid, maxp = maxp, &
        minp = minp, stats = stats, alpha = alpha, controls = controls, &
        settings = settings, info = info, status = status, cov = cov, &
        args = addinfo, err = errmgr, weights = weights)
    if (errmgr%has_error_occurred()) return

    ! End
    return

    ! Memory Error
100 continue
    call report_memory_error("siso_model_fit_least_squares_1", flag, errmgr)
end subroutine

! --------------------
subroutine siso_model_fit_least_squares_2(fcn, x, ic, p, integrator, ind, &
    maxp, minp, stats, alpha, controls, settings, info, status, cov, xc, yc, &
    constraints, weights, args, &
    err)
    !! Attempts to fit a model of a single-intput, single-output (SISO) dynamic 
    !! system by means of an iterative least-squares solver.  The algorithm
    !! computes the solution to the differential equations numerically, and
    !! compares the output to the known solution via a Levenberg-Marquardt
    !! least-squares solver.
    procedure(ode), pointer, intent(in) :: fcn
        !! The routine containing the ODE's being fit.  To communicate 
        !! model parameters and other relevant information, an instance of the
        !! [[model_information]] type is passed to the optional argument of this
        !! routine.  Use the "select type" construct to access this information.
    class(dynamic_system_measurement), intent(in), dimension(:) :: x
        !! An M-element array of arrays with each array containing the measured
        !! input and output of the system being identified.
    real(real64), intent(in), dimension(:,:) :: ic
        !! An M-by-NEQN matrix of initial condition vectors for the NEQN 
        !! equations in fcn, one set for each of the M sets of data in x.
    real(real64), intent(inout), dimension(:) :: p
        !! An N-element array containing an initial guess at the parameters.  
        !! On output, the computed model parameters.
    class(ode_integrator), intent(inout), optional, target :: integrator
        !! The integrator to use when solving the system equations.  If not
        !! supplied, the default integrator will be used.  The default 
        !! integrator is a Runge-Kutta integrator (Dormand-Prince).
    integer(int32), intent(in), optional :: ind
        !! The index of the ODE in fcn providing the output to fit.  If
        !! no value is supplied, a value of 1 will be utilized.
    real(real64), intent(in), optional, dimension(:) :: maxp
        !! An optional N-element array that can be used as upper limits on the 
        !! parameter values. If no upper limit is requested for a particular 
        !! parameter, utilize a very large value. The internal default is to 
        !! utilize huge() as a value.
    real(real64), intent(in), optional, dimension(:) :: minp
        !! An optional N-element array that can be used as lower limits on the 
        !! parameter values. If no lower limit is requested for a particalar 
        !! parameter, utilize a very large magnitude, but negative, value. The 
        !! internal default is to utilize -huge() as a value.
    type(regression_statistics), intent(out), optional, dimension(:) :: stats
        !! An optional N-element array that, if supplied, will be used to 
        !! return statistics about the fit for each parameter.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence 
        !! intervals. The default value is 0.05 such that a 95% confidence 
        !! interval is calculated.
    type(iteration_controls), intent(in), optional :: controls
        !! An optional input providing custom iteration controls.
    type(lm_solver_options), intent(in), optional :: settings
        !! An optional input providing custom settings for the solver.
    type(convergence_info), intent(out), optional :: info
        !! An optional output that can be used to gain information about the 
        !! iterative solution and the nature of the convergence.
    procedure(iteration_update), pointer, intent(in), optional :: status
        !! An optional pointer to a routine that can be used to extract 
        !! iteration information.
    real(real64), intent(out), optional, dimension(:,:) :: cov
        !! An optional N-by-N matrix that, if supplied, will be used to return 
        !! the covariance matrix.
    real(real64), intent(in), optional, dimension(:) :: xc
        !! An optional NC-element array containing the values of the independent 
        !! variable at which the constraint equations are defined.
    real(real64), intent(in), optional, dimension(:) :: yc
        !! An optional NC-element array containing the constraint function 
        !! values at xc.
    procedure(constraint_equations), pointer, optional :: constraints
        !! An optional input, that must be utilized with the xc and yc inputs,
        !! but allows for the implementation of additional constraints on the
        !! solution outside of the differential equations being fitted.  An
        !! example usage would be an additional set of quasi-static tests that
        !! could help identify a stiffness term, for instance.  Other uses of
        !! course can be imagined.
    real(real64), intent(in), optional, dimension(:) :: weights
        !! An optional array containing weighting factors for every equation.
    class(*), intent(inout), optional, target :: args
        !! User-defined information to pass along to fcn.  These arguments,
        !! if supplied, will be passed through to fcn by means of the
        !! [[model_information]] type.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, i1, i2, n, ni, npts, nc, flag
    real(real64), allocatable, dimension(:) :: t, f, ymod, resid
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    type(runge_kutta_45), target :: default_ode_solver
    type(regression_information) :: addinfo
    procedure(regression_function), pointer :: fcnptr
    
    ! Initialization & Error Checking
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(xc) .and. present(yc) .and. present(constraints)) then
        nc = size(xc)
        if (size(yc) /= nc) then
            call report_array_size_error("siso_model_fit_least_squares_2", &
                "yc", nc, size(yc), errmgr)
            return
        end if
    else
        nc = 0
    end if
    n = size(x)
    fcnptr => nlsq_fun

    if (size(ic, 1) /= n) then
        call report_matrix_size_error("siso_model_fit_least_squares_2", "ic", &
            n, size(ic, 2), size(ic, 1), size(ic, 2), errmgr)
        return
    end if

    allocate(addinfo%start_stop(n, 2), stat = flag)
    if (flag /= 0) go to 100
    i1 = 1
    i2 = 0
    npts = 0
    do i = 1, n
        ni = size(x(i)%t)
        npts = npts + ni
        i2 = i2 + ni
        addinfo%start_stop(i, 1) = i1
        addinfo%start_stop(i, 2) = i2
        i1 = i2 + 1
    end do
    npts = npts + nc
    allocate(t(npts), f(npts), addinfo%excitation_data(npts), ymod(npts), &
        resid(npts), stat = flag)
    if (flag /= 0) go to 100

    do i = 1, n
        i1 = addinfo%start_stop(i, 1)
        i2 = addinfo%start_stop(i, 2)
        t(i1:i2) = x(i)%t
        addinfo%excitation_data(i1:i2) = x(i)%input
        f(i1:i2) = x(i)%output
    end do
    if (nc > 0) then
        i1 = i2 + 1
        t(i1:npts) = xc
        f(i1:npts) = yc
    end if

    if (present(integrator)) then
        addinfo%integrator => integrator
    else
        addinfo%integrator => default_ode_solver
    end if

    addinfo%initial_conditions = ic

    if (present(ind)) then
        addinfo%solution_index = ind
    else
        addinfo%solution_index = 1
    end if

    addinfo%ode_routine => fcn

    if (nc > 0) then
        addinfo%uses_constraints = .true.
        addinfo%constraints => constraints
    else
        addinfo%uses_constraints = .false.
    end if

    if (present(args)) addinfo%user_info => args

    ! Solve the system
    call nonlinear_least_squares(fcnptr, t, f, p, ymod, resid, maxp = maxp, &
        minp = minp, stats = stats, alpha = alpha, controls = controls, &
        settings = settings, info = info, status = status, cov = cov, &
        args = addinfo, err = errmgr, weights = weights)
    if (errmgr%has_error_occurred()) return

    ! End
    return

    ! Memory Error
100 continue
    call report_memory_error("siso_model_fit_least_squares_2", flag, errmgr)
end subroutine

! --------------------
subroutine nlsq_fun(t, p, f, check, args)
    !! The routine called by the nonlinear least-squares solver from the
    !! siso_model_fit_least_squares routine.
    real(real64), intent(in), dimension(:) :: t
        !! All of the time values from every data set supplied to the 
        !! least-squares solver.  The args instance contains information on
        !! how to break apart this array.
    real(real64), intent(in), dimension(:) :: p
        !! An array containing the model parameters.
    real(real64), intent(out), dimension(:) :: f
        !! An array, the same size as t, where the function values will be
        !! written.
    logical, intent(out) :: check
        !! Set to true to force a stop to the iteration process; else, set
        !! to false to allow iterations to proceed normally.
    class(*), intent(inout), optional :: args
        !! A regression_information object containing information on how to
        !! arrange and solve the differential equations.

    ! Local Variables
    type(model_information) :: ode_info
    integer(int32) :: i, n, ind, i1, i2
    integer(int32), allocatable, dimension(:,:) :: limits
    real(real64), allocatable, dimension(:) :: excitation, ici
    real(real64), allocatable, dimension(:,:) :: sol, ic
    class(ode_integrator), pointer :: integrator
    type(linear_interpolator), target :: forcing_function
    type(ode_container) :: mdl
    type(errors) :: err
    logical :: uses_constraints
    procedure(constraint_equations), pointer :: constraints
    class(*), pointer :: user_info

    ! Get the supplied information
    select type (args)
    class is (regression_information)
        limits = args%start_stop
        integrator => args%integrator
        ind = args%solution_index + 1   ! +1 accounts for the time vector
        excitation = args%excitation_data
        ic = args%initial_conditions
        mdl%fcn => args%ode_routine
        if (associated(args%user_info)) then
            ode_info%user_info => args%user_info
            user_info => args%user_info
        end if
        uses_constraints = args%uses_constraints
        if (uses_constraints) constraints => args%constraints
    end select

    ! Initialization
    check = .false.
    n = size(limits, 1)
    ode_info%model = p

    ! Cycle over each segment and compute the solution
    do i = 1, n
        ! Get the locations within the array at which the relevant data is
        ! located
        i1 = limits(i, 1)
        i2 = limits(i, 2)

        ! Set up the interpolator for the forcing term
        call forcing_function%initialize(t(i1:i2), excitation(i1:i2), err = err)
        if (err%has_error_occurred()) go to 10
        ode_info%excitation => forcing_function

        ! Solve the ODE's
        if (size(ic, 1) == 1) then
            ici = ic(1,:)
        else
            ici = ic(i,:)
        end if
        call integrator%clear_buffer()
        call integrator%solve(mdl, t(i1:i2), ici, args = ode_info, err = err)
        if (err%has_error_occurred()) go to 10
        sol = integrator%get_solution()
        f(i1:i2) = sol(:,ind)
    end do

    ! Constraints
    if (uses_constraints) then
        if (associated(user_info)) then
            call constraints(t(1:i2), f(1:i2), t(i2+1:), p, f(i2+1:), &
                args = user_info)
        else
            call constraints(t(1:i2), f(1:i2), t(i2+1:), p, f(i2+1:))
        end if
    end if

    ! End
    return

    ! Error Handling
10  continue
    check = .true.  ! terminate iterations
    return
end subroutine

! ------------------------------------------------------------------------------
end module