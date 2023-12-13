module dynamics_error_handling
    use ferror
    use iso_fortran_env
    use diffeq, only : DIFFEQ_INVALID_INPUT_ERROR, &
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
        !! Defins an error associated with a zero-valued frequency.
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

! ------------------------------------------------------------------------------
end module