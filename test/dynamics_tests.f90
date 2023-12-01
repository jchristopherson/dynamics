program main
    use iso_fortran_env, only : int32
    use dynamics_frf_tests
    implicit none

    ! Variables
    logical :: check
    integer(int32) :: flag

    ! Initialization
    flag = 0

    ! Tests
    check = test_fft_based_frf()
    if (.not.check) flag = 1


    ! End
    if (flag /= 0) stop flag
end program