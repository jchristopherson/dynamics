program main
    use iso_fortran_env, only : int32
    use dynamics_frf_tests
    use dynamics_structures_tests
    implicit none

    ! Variables
    logical :: check
    integer(int32) :: flag

    ! Initialization
    flag = 0

    ! Tests
    check = test_frf_sweep()
    if (.not.check) flag = 1

    check = test_proportional_damping_frf()
    if (.not.check) flag = 2

    check = test_modal_response()
    if (.not.check) flag = 3

    check = test_beam2d_shape_functions()
    if (.not.check) flag = 4

    check = test_shape_function_derivatives()
    if (.not.check) flag = 5

    check = test_beam2d_strain_displacement()
    if (.not.check) flag = 6

    check = test_beam2d_stiffness_matrix()
    if (.not.check) flag = 7

    check = test_beam2d_mass_matrix()
    if (.not.check) flag = 8

    check = test_beam2d_ext_force()
    if (.not.check) flag = 9

    check = test_boundary_conditions()
    if (.not.check) flag = 10

    check = test_connectivity_matrix()
    if (.not.check) flag = 11


    ! End
    if (flag /= 0) stop flag
end program