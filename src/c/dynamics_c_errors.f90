! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
module dynamics_c_errors
    ! Non-fatal error reporting for the C API.  Argument errors detected by the
    ! C wrappers are recorded per thread and forwarded to an optional
    ! user-supplied handler rather than terminating the process.
    use iso_c_binding
    use iso_fortran_env
    implicit none
    private
    public :: c_api_error
    public :: c_error_handler
    public :: c_set_error_handler
    public :: c_get_last_error
    public :: c_get_last_error_message
    public :: c_clear_error

    interface
        subroutine c_error_handler(code, message, user_data) &
            bind(C, name = "c_error_handler")
            use iso_c_binding
            integer(c_int), intent(in), value :: code
            character(kind = c_char), intent(in) :: message(*)
            type(c_ptr), intent(in), value :: user_data
        end subroutine
    end interface

    integer(c_int), parameter :: MESSAGE_LENGTH = 256

    integer(c_int) :: last_error = 0
    character(len = MESSAGE_LENGTH) :: last_message = ""
    type(c_funptr) :: handler = c_null_funptr
    type(c_ptr) :: handler_data = c_null_ptr
    !$omp threadprivate(last_error, last_message, handler, handler_data)

contains

function c_api_error(condition, code, message) result(rst)
    ! Records and reports an error if condition is true.  Returns condition so
    ! the caller can return early.
    logical, intent(in) :: condition
    integer(int32), intent(in) :: code
    character(len = *), intent(in) :: message
    logical :: rst

    procedure(c_error_handler), pointer :: fcn
    character(kind = c_char, len = MESSAGE_LENGTH + 1) :: cmsg

    rst = condition
    if (.not.condition) return
    last_error = int(code, c_int)
    last_message = message
    if (c_associated(handler)) then
        cmsg = trim(last_message) // c_null_char
        call c_f_procpointer(handler, fcn)
        call fcn(last_error, cmsg, handler_data)
    end if
end function

! ------------------------------------------------------------------------------
subroutine c_set_error_handler(fcn, user_data) &
    bind(C, name = "c_set_error_handler")
    type(c_funptr), intent(in), value :: fcn
    type(c_ptr), intent(in), value :: user_data
    handler = fcn
    handler_data = user_data
end subroutine

! ------------------------------------------------------------------------------
function c_get_last_error() result(rst) bind(C, name = "c_get_last_error")
    integer(c_int) :: rst
    rst = last_error
end function

! ------------------------------------------------------------------------------
function c_get_last_error_message(n, buffer) result(rst) &
    bind(C, name = "c_get_last_error_message")
    integer(c_int), intent(in), value :: n
    character(kind = c_char), intent(out) :: buffer(*)
    integer(c_int) :: rst

    integer(c_int) :: i, m

    rst = int(len_trim(last_message), c_int)
    if (n < 1) return
    m = min(rst, n - 1)
    do i = 1, m
        buffer(i) = last_message(i:i)
    end do
    buffer(m + 1) = c_null_char
end function

! ------------------------------------------------------------------------------
subroutine c_clear_error() bind(C, name = "c_clear_error")
    last_error = 0
    last_message = ""
end subroutine

! ------------------------------------------------------------------------------
end module
