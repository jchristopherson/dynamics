module dynamics_solver_reporting
    !! Provides mechanisms for reporting solver status. Reporting is optional:
    !! every progress callback is an optional argument, and solvers do no
    !! reporting work when it is omitted. Static solvers report accepted
    !! nonlinear iterates, while time-history solvers report completed steps.
    use iso_fortran_env, only : real64
    implicit none
    private
    public :: static_solver_progress
    public :: time_history_solver_progress
    public :: static_solver_callback
    public :: time_history_solver_callback

    type :: static_solver_progress
        !! Describes an accepted iterate of a nonlinear static solver.
        integer :: iteration = 0
            !! The number of accepted iterations completed so far.
        real(real64) :: residual_norm = 0.0d0
            !! The Euclidean norm of the free residual at the accepted iterate.
        real(real64) :: tolerance = 0.0d0
            !! The residual norm required for convergence.
        real(real64) :: step_fraction = 1.0d0
            !! The fraction of the Newton correction accepted after
            !! backtracking. A value of one indicates a full Newton step.
    end type

    type :: time_history_solver_progress
        !! Describes a completed time step of a time-history solver.
        integer :: step = 0
            !! The number of time steps completed in the current solve.
        integer :: total_steps = 0
            !! The number of time steps requested by the current solve.
        real(real64) :: time = 0.0d0
            !! The elapsed time since the start of the current solve.
        real(real64), allocatable, dimension(:) :: displacement
            !! A copy of the displacement at the completed time point.
        real(real64), allocatable, dimension(:) :: velocity
            !! A copy of the velocity at the completed time point.
        real(real64), allocatable, dimension(:) :: acceleration
            !! A copy of the acceleration at the completed time point.
    end type

    abstract interface
        subroutine static_solver_callback(progress, args)
            !! Receives progress from a nonlinear static solver.
            import :: static_solver_progress
            type(static_solver_progress), intent(in) :: progress
                !! The accepted iterate being reported.
            class(*), intent(in), optional :: args
                !! Optional user-supplied data.
        end subroutine

        subroutine time_history_solver_callback(progress, args)
            !! Receives progress from a time-history solver.
            import :: time_history_solver_progress
            type(time_history_solver_progress), intent(in) :: progress
                !! The completed time step being reported.
            class(*), intent(in), optional :: args
                !! Optional user-supplied data.
        end subroutine
    end interface
end module