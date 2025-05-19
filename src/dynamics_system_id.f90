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
        real(real64), allocatable, dimension(:) :: initial_conditions
            !! The initial conditions array to pass to the ODE solver.
        real(real64), allocatable, dimension(:) :: excitation_data
            !! An array of excitation data to pass to the ODE solver.
        procedure(ode), pointer, nopass :: ode_routine
            !! The routine, defined by the calling code, containing the ODE's
            !! to solve.
    end type

contains
! ------------------------------------------------------------------------------
subroutine siso_model_fit_least_squares(x, p, integrator, err)
    !! Attempts to fit a model of a single-intput, single-output (SISO) dynamic 
    !! system by means of an iterative least-squares solver
    class(dynamic_system_measurement), intent(in), dimension(:) :: x
        !! An M-element array of arrays with each array containing the measured
        !! input and output of the system being identified.
    real(real64), intent(inout), dimension(:) :: p
        !! An array containing an initial guess at the parameters.  On output,
        !! the computed model parameters.
    class(ode_integrator), intent(inout), optional, target :: integrator
        !! The integrator to use when solving the system equations.
    class(errors), intent(inout), optional, target :: err
        !! An error handling object.

    ! Local Variables
    integer(int32) :: i, n
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    class(ode_integrator), pointer :: solver_ode
    type(rosenbrock), target :: default_ode_solver
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(integrator)) then
        solver_ode => integrator
    else
        solver_ode => default_ode_solver
    end if
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
    real(real64), allocatable, dimension(:) :: excitation, ic
    real(real64), allocatable, dimension(:,:) :: sol
    class(ode_integrator), pointer :: integrator
    type(linear_interpolator), target :: forcing_function
    type(ode_container) :: mdl
    type(errors) :: err

    ! Get the supplied information
    select type (args)
    class is (regression_information)
        limits = args%start_stop
        integrator => args%integrator
        ind = args%solution_index + 1   ! +1 accounts for the time vector
        excitation = args%excitation_data
        ic = args%initial_conditions
        mdl%fcn => args%ode_routine
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
        call integrator%clear_buffer()
        call integrator%solve(mdl, t(i1:i2), ic, args = ode_info, err = err)
        if (err%has_error_occurred()) go to 10
        sol = integrator%get_solution()
        f(i1:i2) = sol(:,ind)
    end do

    ! End
    return

    ! Error Handling
10  continue
    check = .true.  ! terminate iterations
    return
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module