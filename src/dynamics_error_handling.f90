module dynamics_error_handling
    use ferror
    use iso_fortran_env
    use diffeq_errors, only : DIFFEQ_INVALID_INPUT_ERROR, &
        DIFFEQ_MEMORY_ALLOCATION_ERROR, DIFFEQ_NULL_POINTER_ERROR
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
contains
! ------------------------------------------------------------------------------
    subroutine report_null_forcing_routine_error(name, err)
        !! Reports a null forcing routine pointer error.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.
        
        ! Report the error
        call err%report_error(name, &
            "No forcing function routine was supplied.", &
            DYN_NULL_POINTER_ERROR)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_memory_error(name, flag, err)
        !! Reports a memory allocation error.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        integer(int32), intent(in) :: flag
            !! The flag returned from the allocate statement.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) "Memory allocation error flag ", flag, "."
        call err%report_error(name, trim(errmsg), DYN_MEMORY_ERROR)

        ! Formatting
    100 format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_nonsquare_mass_matrix_error(name, m, n, err)
        !! Reports an error relating to a non-square mass matrix.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        integer(int32), intent(in) :: m
            !! The number of rows found in the mass matrix.
        integer(int32), intent(in) :: n
            !! The number of columns found in the mass matrix.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg
        
        ! Report the error
        write(errmsg, 100) "The mass matrix is not square.  " // &
            "It was found to be ", m, "-by-", n, "."
        call err%report_error(name, trim(errmsg), DYN_MATRIX_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_nonsquare_stiffness_matrix_error(name, m, n, err)
        !! Reports an error relating to a non-square stiffness matrix.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        integer(int32), intent(in) :: m
            !! The number of rows found in the stiffness matrix.
        integer(int32), intent(in) :: n
            !! The number of columns found in the stiffness matrix.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) "The stiffness matrix is not square.  " // &
            "It was found to be ", m, "-by-", n, "."
        call err%report_error(name, trim(errmsg), DYN_MATRIX_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_nonsquare_matrix_error(name, var, m, n, err)
        !! Reports an error relating to a non-square matrix.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        character(len = *), intent(in) :: var
            !! The name of the offending variable.
        integer(int32), intent(in) :: m
            !! The number of rows found in the matrix.
        integer(int32), intent(in) :: n
            !! The number of columns found in the matrix.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) "Matrix " // var // " is not square.  " // &
            "It was found to be ", m, "-by-", n, "."
        call err%report_error(name, trim(errmsg), DYN_MATRIX_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_matrix_size_mismatch_error(name, mtx1, mtx2, m1, n1, &
        m2, n2, err)
        !! Reports a mismatch in matrix sizes.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        character(len = *), intent(in) :: mtx1
            !! The name of the first matrix.
        character(len = *), intent(in) :: mtx2
            !! The name of the second matrix.
        integer(int32), intent(in) :: m1
            !! The number of rows in the first matrix.
        integer(int32), intent(in) :: n1
            !! The number of columns in the first matrix.
        integer(int32), intent(in) :: m2
            !! The number of rows in the second matrix.
        integer(int32), intent(in) :: n2
            !! The number of columns in the second matrix.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg
        
        ! Report the error
        write(errmsg, 100) "The size of the " // mtx1 // " matrix (", m1, &
            "-by-", n1, ") does not match the size of the " // mtx2, &
            " matrix (", m2, "-by-", n2, ")."
        call err%report_error(name, trim(errmsg), DYN_MATRIX_SIZE_ERROR)

        ! Formatting
    100 format(A, I0, A, I0, A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_zero_valued_frequency_error(name, index, err)
        !! Reports an error associated with a zero-valued frequency value.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        integer(int32), intent(in) :: index
            !! The array index at which the zero-valued frequency was found.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) "A zero-valued frequency was found at index ", &
            index, "."
        call err%report_error(name, trim(errmsg), &
            DYN_ZERO_VALUED_FREQUENCY_ERROR)

        ! Formatting
    100 format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_generic_counting_error(name, str1, val, str2, flag, err)
        !! A generic error reporting routine.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        character(len = *), intent(in) :: str1
            !! The first string.
        integer(int32), intent(in) :: val
            !! The integer value.
        character(len = *), intent(in) :: str2
            !! The second string.
        integer(int32), intent(in) :: flag
            !! The error flag.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 512) :: errmsg

        ! Report the error
        write(errmsg, 100) str1, val, str2
        call err%report_error(name, trim(errmsg), flag)

        ! Formatting
    100 format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_zero_difference_error(name, var1, val1, var2, val2, &
        flag, err)
        !! Reports a zero-difference between two variables where a non-zero
        !! difference was expected.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        character(len = *), intent(in) :: var1
            !! The name of the first variable.
        real(real64), intent(in) :: val1
            !! The value of the first variable.
        character(len = *), intent(in) :: var2
            !! The name of the second variable.
        real(real64), intent(in) :: val2
            !! The value of the second variable.
        integer(int32), intent(in) :: flag
            !! The error flag.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) "A non-zero difference between " // var1 // &
            " (", val1, "), and " // var2 // " (", val2, ") was expected."
        call err%report_error(name, trim(errmsg), flag)

        ! Formatting
    100 format(A, F0.4, A, F0.4, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_overconstraint_error(name, err)
        !! Reports an overconstraint error.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Report the error
        call err%report_error(name, "The model is overconstrained.", &
            DYN_CONSTRAINT_ERROR)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_array_index_out_of_bounds_error(name, var, ind, sz, err)
        !! Reports an array index-out-of-bounds error.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        character(len = *), intent(in) :: var
            !! The name of the offending variable.
        integer(int32), intent(in) :: ind
            !! The offending index.
        integer(int32), intent(in) :: sz
            !! The array size.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) &
            "Index ", ind, " is outside the bounds of array " // var // &
            ", which is of size ", sz, "."
        call err%report_error(name, trim(errmsg), DYN_INDEX_OUT_OF_RANGE)

        ! Formatting
    100 format(A, I0, A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine report_nonmonotonic_array_error(name, var, ind, err)
        !! Reports a nonmonotonic array error.
        character(len = *), intent(in) :: name
            !! The name of the routine in which the error was found.
        character(len = *), intent(in) :: var
            !! The name of the offending variable.
        integer(int32), intent(in) :: ind
            !! The index of the occurrence of nonmonotonicity.
        class(errors), intent(inout) :: err
            !! An errors-based object that if provided can be used to retrieve 
            !! information relating to any errors encountered during execution.

        ! Local Variables
        character(len = 256) :: errmsg

        ! Report the error
        write(errmsg, 100) "Array " // var // &
            " was found to be nonmonotonic at index ", ind, "."
        call err%report_error(name, trim(errmsg), DYN_NONMONOTONIC_ARRAY_ERROR)

        ! Formatting
    100 format(A, I0, A)
    end subroutine

! ------------------------------------------------------------------------------
end module