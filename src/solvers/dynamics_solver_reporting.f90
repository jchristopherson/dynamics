module dynamics_solver_reporting
    !! Provides mechanisms for reporting solver status.
    use iso_fortran_env, only : real64
    implicit none
    private

    ! TO DO: Define subroutine interfaces for user callbacks that can be
    ! used to gain information about solver progress.  There should be callbacks
    ! for both static and time history solvers as their requirements are 
    ! different.  It may be worth having a few custom types that contain
    ! solution information to keep the callback signatures manageable and
    ! extensible
end module