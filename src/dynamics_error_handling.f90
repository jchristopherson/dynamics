module dynamics_error_handling
    use iso_fortran_env
    use diffeq_errors, only : DIFFEQ_INVALID_INPUT_ERROR, &
        DIFFEQ_MEMORY_ALLOCATION_ERROR, DIFFEQ_NULL_POINTER_ERROR
    use fstats_errors, only : FS_UNDERDEFINED_PROBLEM_ERROR, &
        FS_TOLERANCE_TOO_SMALL_ERROR, FS_TOO_FEW_ITERATION_ERROR
    use nonlin_error_handling, only : NL_CONVERGENCE_ERROR
    implicit none

    integer(int32), parameter :: DYN_MEMORY_ERROR = DIFFEQ_MEMORY_ALLOCATION_ERROR
        !! Defines an error associated with memory allocations.
    integer(int32), parameter :: DYN_NULL_POINTER_ERROR = DIFFEQ_NULL_POINTER_ERROR
        !! Defines an error associated with a null pointer.
    integer(int32), parameter :: DYN_INVALID_INPUT_ERROR = DIFFEQ_INVALID_INPUT_ERROR
        !! Defines an error associated with an invalid input.
    integer(int32), parameter :: DYN_MATRIX_SIZE_ERROR = 100100
        !! Defines an error associated with an incorrectly sized matrix.
    integer(int32), parameter :: DYN_ZERO_VALUED_FREQUENCY_ERROR = 100101
        !! Defines an error associated with a zero-valued frequency.
    integer(int32), parameter :: DYN_CONSTRAINT_ERROR = 100102
        !! Defines a constraint-related error.
    integer(int32), parameter :: DYN_INDEX_OUT_OF_RANGE = 100103
        !! Defines an index out of range error.
    integer(int32), parameter :: DYN_NONMONOTONIC_ARRAY_ERROR = 100104
        !! Defines an error related to an array being nonmonotonic.
    integer(int32), parameter :: DYN_ARRAY_SIZE_ERROR = 100105
        !! Defines an error for an improperly sized array.
    integer(int32), parameter :: DYN_UNDERDEFINED_PROBLEM_EROR = FS_UNDERDEFINED_PROBLEM_ERROR
        !! Defines an error for an underdefined problem.
    integer(int32), parameter :: DYN_TOLERANCE_TOO_SMALL_ERROR = FS_TOLERANCE_TOO_SMALL_ERROR
        !! Defines an error related to the request of a too small tolerance 
        !! value.
    integer(int32), parameter :: DYN_TOO_FEW_ITERATIONS_ERROR = FS_TOO_FEW_ITERATION_ERROR
        !! Defines an error when too few iterations were allowed.
    integer(int32), parameter :: DYN_CONVERGENCE_ERROR = NL_CONVERGENCE_ERROR
        !! Defines an error related to convergence issues.

end module