module dynamics_rotation_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_rotation
    implicit none

contains
! ------------------------------------------------------------------------------
function test_quaternion_init_1() result(rst)
    ! Variables
    logical :: rst
    type(quaternion) :: q
    real(real64) :: x(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)

    ! Test
    if (.not.assert(x(1), q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 1"
    end if

    if (.not.assert(x(2), q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 2"
    end if

    if (.not.assert(x(3), q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 3"
    end if

    if (.not.assert(x(4), q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_init_2() result(rst)
    ! Variables
    logical :: rst
    type(quaternion) :: q
    real(real64) :: angle, axis(3), w, x, y, z

    ! Initialization
    rst = .true.
    call random_number(angle)
    call random_number(axis)
    w = cos(0.5d0 * angle)
    x = sin(0.5d0 * angle) * axis(1)
    y = sin(0.5d0 * angle) * axis(2)
    z = sin(0.5d0 * angle) * axis(3)
    q = quaternion(angle, axis)

    ! Test
    if (.not.assert(w, q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 1"
    end if

    if (.not.assert(x, q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 2"
    end if

    if (.not.assert(y, q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 3"
    end if

    if (.not.assert(z, q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 4"
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module