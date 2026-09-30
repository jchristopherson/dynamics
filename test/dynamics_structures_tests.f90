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
module dynamics_structures_tests
    use iso_fortran_env
    use dynamics
    use fortran_test_helper
    implicit none

contains
! ------------------------------------------------------------------------------
    pure function beam2d_n1(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.5d0 * (1.0d0 - s)
    end function
    
    pure function beam2d_n2(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.25d0 * (1.0d0 - s)**2 * (2.0d0 + s)
    end function
    
    pure function beam2d_n3(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.25d0 * (1.0d0 - s)**2 * (1.0d0 + s)
    end function
    
    pure function beam2d_n4(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.5d0 * (1.0d0 + s)
    end function
    
    pure function beam2d_n5(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.25d0 * (2.0d0 - s) * (1.0d0 + s)**2
    end function
    
    pure function beam2d_n6(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.25d0 * (1.0d0 + s)**2 * (s - 1.0d0)
    end function
    
    pure function beam2d_shape_fcn_mtx(s, l) result(rst)
        real(real64), intent(in) :: s, l
        real(real64) :: rst(2, 6)
    
        real(real64), parameter :: z = 0.0d0
        real(real64) :: n1, n2, n3, n4, n5, n6
    
        n1 = beam2d_n1(s)
        n2 = beam2d_n2(s)
        n3 = 0.5d0 * l * beam2d_n3(s)
        n4 = beam2d_n4(s)
        n5 = beam2d_n5(s)
        n6 = 0.5d0 * l * beam2d_n6(s)
    
        rst = reshape([n1, z, z, n2, z, n3, n4, z, z, n5, z, n6], [2, 6])
    end function
    
    pure function beam2d_strain_disp_matrix(s, l) result(rst)
        real(real64), intent(in) :: s, l
        real(real64) :: rst(2, 6)
    
        rst = reshape([ &
            -1.0d0 / l, 0.0d0, &
            0.0d0, 6.0d0 * s / (l**2), &
            0.0d0, (3.0d0 * s - 1.0d0) / l, &
            1.0d0 / l, 0.0d0, &
            0.0d0, -6.0d0 * s / (l**2), &
            0.0d0, (3.0d0 * s + 1.0d0) / l &
        ], [2, 6])
    end function
    
    ! ------------------------------------------------------------------------------
    function test_beam2d_shape_functions() result(rst)
        logical :: rst
    
        ! Parameters
        real(real64), parameter :: s1(1) = [-sqrt(3.0d0) / 3.0d0]
        real(real64), parameter :: s2(1) = [0.0d0]
        real(real64), parameter :: s3(1) = [-s1]
    
        ! Local Variables
        real(real64) :: l, x1, y1, x2, y2
        type(beam_element_2d) :: e
    
        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(y1)
        call random_number(y2)
        e%node_1%index = 1
        e%node_1%x = x1
        e%node_1%y = y1
        e%node_1%z = 0.0d0
        e%node_2%index = 2
        e%node_2%x = x2
        e%node_2%y = y2
        e%node_2%z = 0.0d0
        e%node_1%dof = 3
        e%node_2%dof = 3
        l = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
        ! Tests - #1
        if (.not.assert( &
            e%evaluate_shape_function(1, s1), &
            beam2d_n1(s1(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -1"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(1, s2), &
            beam2d_n1(s2(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -2"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(1, s3), &
            beam2d_n1(s3(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -3"
        end if
    
        ! Tests - #2
        if (.not.assert( &
            e%evaluate_shape_function(2, s1), &
            beam2d_n2(s1(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -4"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(2, s2), &
            beam2d_n2(s2(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -5"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(2, s3), &
            beam2d_n2(s3(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -6"
        end if
    
        ! Tests - #3
        if (.not.assert( &
            e%evaluate_shape_function(3, s1), &
            beam2d_n3(s1(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -7"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(3, s2), &
            beam2d_n3(s2(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -8"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(3, s3), &
            beam2d_n3(s3(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -9"
        end if
    
        ! Tests - #4
        if (.not.assert( &
            e%evaluate_shape_function(4, s1), &
            beam2d_n4(s1(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -10"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(4, s2), &
            beam2d_n4(s2(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -11"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(4, s3), &
            beam2d_n4(s3(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -12"
        end if
    
        ! Tests - #5
        if (.not.assert( &
            e%evaluate_shape_function(5, s1), &
            beam2d_n5(s1(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -13"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(5, s2), &
            beam2d_n5(s2(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -14"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(5, s3), &
            beam2d_n5(s3(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -15"
        end if
    
        ! Tests - #6
        if (.not.assert( &
            e%evaluate_shape_function(6, s1), &
            beam2d_n6(s1(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -16"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(6, s2), &
            beam2d_n6(s2(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -17"
        end if
    
        if (.not.assert( &
            e%evaluate_shape_function(6, s3), &
            beam2d_n6(s3(1)))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -18"
        end if
    
        ! Shape function matrix tests
        if (.not.assert( &
            e%shape_function_matrix(s1), &
            beam2d_shape_fcn_mtx(s1(1), l))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -19"
        end if
        if (.not.assert( &
            e%shape_function_matrix(s2), &
            beam2d_shape_fcn_mtx(s2(1), l))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -20"
        end if
        if (.not.assert( &
            e%shape_function_matrix(s3), &
            beam2d_shape_fcn_mtx(s3(1), l))) &
        then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_shape_functions -21"
        end if
    end function

! ------------------------------------------------------------------------------
    function dN5ds(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = 0.75d0 * (1.0d0 - s**2)
    end function

    function dN5ds2(s) result(rst)
        real(real64), intent(in) :: s
        real(real64) :: rst
        rst = -1.5d0 * s
    end function

    function test_shape_function_derivatives() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6
        real(real64), parameter :: s1(1) = [-sqrt(3.0d0) / 3.0d0]
        real(real64), parameter :: s2(1) = [0.0d0]
        real(real64), parameter :: s3(1) = [-s1]

        ! Local Variables
        real(real64) :: x1, y1, x2, y2, deriv, ans
        type(beam_element_2d) :: e
    
        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(y1)
        call random_number(y2)
        e%node_1%index = 1
        e%node_1%x = x1
        e%node_1%y = y1
        e%node_1%z = 0.0d0
        e%node_2%index = 2
        e%node_2%x = x2
        e%node_2%y = y2
        e%node_2%z = 0.0d0
        e%node_1%dof = 3
        e%node_2%dof = 3

        ! Test the first derivatives
        deriv = shape_function_derivative(5, e, s1, 1)
        ans = dN5ds(s1(1))
        if (.not.assert(deriv, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_shape_function_derivatives -1"
        end if

        deriv = shape_function_derivative(5, e, s2, 1)
        ans = dN5ds(s2(1))
        if (.not.assert(deriv, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_shape_function_derivatives -2"
        end if

        deriv = shape_function_derivative(5, e, s3, 1)
        ans = dN5ds(s3(1))
        if (.not.assert(deriv, ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_shape_function_derivatives -3"
        end if

        ! Test the second derivatives
        deriv = shape_function_second_derivative(5, e, s1, 1)
        ans = dN5ds2(s1(1))
        if (.not.assert(deriv, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_shape_function_derivatives -4"
        end if

        deriv = shape_function_second_derivative(5, e, s2, 1)
        ans = dN5ds2(s2(1))
        if (.not.assert(deriv, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_shape_function_derivatives -5"
        end if

        deriv = shape_function_second_derivative(5, e, s3, 1)
        ans = dN5ds2(s3(1))
        if (.not.assert(deriv, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_shape_function_derivatives -6"
        end if
    end function
    
! ------------------------------------------------------------------------------
    function test_beam2d_strain_displacement() result(rst)
        ! Arguments
        logical :: rst
    
        ! Parameters
        real(real64), parameter :: tol = 5.0d-6
        real(real64), parameter :: s1(1) = [-sqrt(3.0d0) / 3.0d0]
        real(real64), parameter :: s2(1) = [0.0d0]
        real(real64), parameter :: s3(1) = [-s1]
    
        ! Local Variables
        real(real64) :: l, x1, y1, x2, y2
        real(real64), allocatable, dimension(:,:) :: b, ans
        type(beam_element_2d) :: e
    
        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(y1)
        call random_number(y2)
        e%node_1%index = 1
        e%node_1%x = x1
        e%node_1%y = y1
        e%node_1%z = 0.0d0
        e%node_2%index = 2
        e%node_2%x = x2
        e%node_2%y = y2
        e%node_2%z = 0.0d0
        e%node_1%dof = 3
        e%node_2%dof = 3
        l = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
        ! Tests
        b = e%strain_displacement_matrix(s1)
        ans = beam2d_strain_disp_matrix(s1(1), l)
        if (.not.assert(b, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_strain_displacement -1"
        end if
    
        b = e%strain_displacement_matrix(s2)
        ans = beam2d_strain_disp_matrix(s2(1), l)
        if (.not.assert(b, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_strain_displacement -2"
        end if
    
        b = e%strain_displacement_matrix(s3)
        ans = beam2d_strain_disp_matrix(s3(1), l)
        if (.not.assert(b, ans, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_strain_displacement -3"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_stress() result(rst)
        logical :: rst

        real(real64), parameter :: tol = 1.0d-10
        real(real64), parameter :: area = 2.5d0
        real(real64), parameter :: modulus = 3.0d7
        real(real64), parameter :: displacement = 2.0d-3
        real(real64), parameter :: length = 5.0d0
        real(real64), parameter :: cosine = 3.0d0 / 5.0d0
        real(real64), parameter :: sine = 4.0d0 / 5.0d0
        real(real64), parameter :: s(1) = [0.25d0]

        real(real64) :: u(6), expected(2)
        real(real64), allocatable, dimension(:) :: stress
        type(beam_element_2d) :: e

        rst = .true.
        e%node_1 = node(1, 3, 1.0d0, 2.0d0, 0.0d0)
        e%node_2 = node(2, 3, 4.0d0, 6.0d0, 0.0d0)
        e%area = area
        e%moment_of_inertia = 0.75d0
        e%material%modulus = modulus

        u = [0.0d0, 0.0d0, 0.0d0, displacement * cosine, &
            displacement * sine, 0.0d0]
        expected = [area * modulus * displacement / length, 0.0d0]
        stress = e%stress(u, s)

        if (.not.assert(stress, expected, tol * maxval(abs(expected)))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_stress -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_strain() result(rst)
        logical :: rst

        real(real64), parameter :: tol = 1.0d-12
        real(real64), parameter :: displacement = 2.0d-3
        real(real64), parameter :: length = 5.0d0
        real(real64), parameter :: cosine = 3.0d0 / 5.0d0
        real(real64), parameter :: sine = 4.0d0 / 5.0d0
        real(real64), parameter :: s(1) = [0.25d0]

        real(real64) :: u(6), expected(2)
        real(real64), allocatable, dimension(:) :: strain
        type(beam_element_2d) :: e

        rst = .true.
        e%node_1 = node(1, 3, 1.0d0, 2.0d0, 0.0d0)
        e%node_2 = node(2, 3, 4.0d0, 6.0d0, 0.0d0)

        u = [0.0d0, 0.0d0, 0.0d0, displacement * cosine, &
            displacement * sine, 0.0d0]
        expected = [displacement / length, 0.0d0]
        strain = e%strain(u, s)

        if (.not.assert(strain, expected, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_strain -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_bending_stress() result(rst)
        logical :: rst

        real(real64), parameter :: length = 5.0d0
        real(real64), parameter :: area = 2.5d0
        real(real64), parameter :: moi = 0.75d0
        real(real64), parameter :: modulus = 3.0d7
        real(real64), parameter :: curvature = 2.0d-3
        real(real64), parameter :: tol = 1.0d-12
        real(real64), parameter :: s(1) = [0.25d0]

        real(real64) :: displacement(6), expected(2)
        real(real64), allocatable :: stress(:)
        type(material) :: mat
        type(beam_element_2d) :: e

        rst = .true.
        mat = material(modulus, 0.3d0, 1.0d0)
        e = beam_element_2d(mat, area, moi, &
            node(1, 3, 0.0d0, 0.0d0, 0.0d0), &
            node(2, 3, length, 0.0d0, 0.0d0))
        displacement = [0.0d0, 0.0d0, 0.0d0, &
            0.0d0, 0.5d0 * curvature * length**2, curvature * length]
        expected = [0.0d0, modulus * moi * curvature]
        stress = e%stress(displacement, s)

        if (.not.assert(stress, expected, tol * maxval(abs(expected)))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_bending_stress -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_bending_strain() result(rst)
        logical :: rst

        real(real64), parameter :: length = 5.0d0
        real(real64), parameter :: curvature = 2.0d-3
        real(real64), parameter :: s(1) = [0.25d0]

        real(real64) :: displacement(6), expected(2)
        real(real64), allocatable :: strain(:)
        type(material) :: mat
        type(beam_element_2d) :: e

        rst = .true.
        mat = material(3.0d7, 0.3d0, 1.0d0)
        e = beam_element_2d(mat, 2.5d0, 0.75d0, &
            node(1, 3, 0.0d0, 0.0d0, 0.0d0), &
            node(2, 3, length, 0.0d0, 0.0d0))
        displacement = [0.0d0, 0.0d0, 0.0d0, &
            0.0d0, 0.5d0 * curvature * length**2, curvature * length]
        expected = [0.0d0, curvature]
        strain = e%strain(displacement, s)

        if (.not.assert(strain, expected)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_bending_strain -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_internal_results() result(rst)
        logical :: rst

        real(real64), parameter :: length = 5.0d0
        real(real64), parameter :: moi = 0.75d0
        real(real64), parameter :: modulus = 3.0d7
        real(real64), parameter :: coefficient = 1.0d-4
        real(real64), parameter :: s(1) = [0.25d0]
        real(real64), parameter :: x = 0.5d0 * length * (s(1) + 1.0d0)
        real(real64), parameter :: tol = 1.0d-12

        real(real64) :: displacement(6), expected_moment, expected_shear
        type(material) :: mat
        type(beam_element_2d) :: e

        rst = .true.
        mat = material(modulus, 0.3d0, 1.0d0)
        e = beam_element_2d(mat, 2.5d0, moi, &
            node(1, 3, 0.0d0, 0.0d0, 0.0d0), &
            node(2, 3, length, 0.0d0, 0.0d0))
        displacement = [0.0d0, 0.0d0, 0.0d0, &
            0.0d0, coefficient * length**3, 3.0d0 * coefficient * length**2]
        expected_moment = modulus * moi * 6.0d0 * coefficient * x
        expected_shear = modulus * moi * 6.0d0 * coefficient

        if (.not.assert(e%bending_moment(displacement, s), expected_moment, &
                tol * abs(expected_moment))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_internal_results -1"
        end if
        if (.not.assert(e%shear_force(displacement, s), expected_shear, &
                tol * abs(expected_shear))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_internal_results -2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_nodally_averaged_stress() result(rst)
        logical :: rst

        real(real64), parameter :: tol = 1.0d-12
        real(real64) :: displacement(9), expected(2,3)
        real(real64), allocatable, dimension(:,:) :: stress
        type(beam_element_2d) :: elements(2)
        type(material) :: mat
        type(node) :: nodes(3)

        rst = .true.
        mat%modulus = 10.0d0
        mat%poissons_ratio = 0.3d0
        mat%density = 1.0d0
        nodes = [ &
            node(10, 3, 0.0d0, 0.0d0, 0.0d0), &
            node(20, 3, 1.0d0, 0.0d0, 0.0d0), &
            node(30, 3, 2.0d0, 0.0d0, 0.0d0) &
        ]
        elements(1) = beam_element_2d(mat, 2.0d0, 1.0d0, &
            nodes(1), nodes(2))
        elements(2) = beam_element_2d(mat, 2.0d0, 1.0d0, &
            nodes(2), nodes(3))
        displacement = [ &
            0.0d0, 0.0d0, 0.0d0, &
            0.1d0, 0.0d0, 0.0d0, &
            0.3d0, 0.0d0, 0.0d0 &
        ]
        expected = reshape([ &
            2.0d0, 0.0d0, &
            3.0d0, 0.0d0, &
            4.0d0, 0.0d0 &
        ], [2, 3])

        stress = nodally_averaged_stress(elements, nodes, displacement)

        if (.not.assert(stress, expected, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_nodally_averaged_stress -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_nodally_averaged_strain() result(rst)
        logical :: rst

        real(real64), parameter :: tol = 1.0d-12
        real(real64) :: displacement(9), expected(2,3)
        real(real64), allocatable, dimension(:,:) :: strain
        type(beam_element_2d) :: elements(2)
        type(material) :: mat
        type(node) :: nodes(3)

        rst = .true.
        mat%modulus = 10.0d0
        mat%poissons_ratio = 0.3d0
        mat%density = 1.0d0
        nodes = [ &
            node(10, 3, 0.0d0, 0.0d0, 0.0d0), &
            node(20, 3, 1.0d0, 0.0d0, 0.0d0), &
            node(30, 3, 2.0d0, 0.0d0, 0.0d0) &
        ]
        elements(1) = beam_element_2d(mat, 2.0d0, 1.0d0, &
            nodes(1), nodes(2))
        elements(2) = beam_element_2d(mat, 2.0d0, 1.0d0, &
            nodes(2), nodes(3))
        displacement = [ &
            0.0d0, 0.0d0, 0.0d0, &
            0.1d0, 0.0d0, 0.0d0, &
            0.3d0, 0.0d0, 0.0d0 &
        ]
        expected = reshape([ &
            0.1d0, 0.0d0, &
            0.15d0, 0.0d0, &
            0.2d0, 0.0d0 &
        ], [2, 3])

        strain = nodally_averaged_strain(elements, nodes, displacement)

        if (.not.assert(strain, expected, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_nodally_averaged_strain -1"
        end if
    end function
    
    ! ------------------------------------------------------------------------------
    function test_beam2d_stiffness_matrix() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6
    
        ! Local Variables
        real(real64) :: l, x1, y1, x2, y2, w, h
        real(real64), allocatable, dimension(:,:) :: k, ans, T
        type(beam_element_2d) :: e
    
        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(y1)
        call random_number(y2)
        e%node_1%index = 1
        e%node_1%x = x1
        e%node_1%y = y1
        e%node_1%z = 0.0d0
        e%node_2%index = 2
        e%node_2%x = x2
        e%node_2%y = y2
        e%node_2%z = 0.0d0
        e%node_1%dof = 3
        e%node_2%dof = 3
        l = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        T = e%rotation_matrix()
    
        ! Define the cross-sectional properties
        call random_number(w)
        call random_number(h)
        e%area = w * h
        e%moment_of_inertia = w * h**3 / 12.0d0
    
        ! Define the material properties
        e%material%density = 0.101d0 / 3.86d2
        e%material%modulus = 10.0d6
        e%material%poissons_ratio = 0.33d0
    
        ! Define the solution
        allocate(ans(6, 6), source = 0.0d0)
        ans(1,1) = e%area * e%material%modulus / l
        ans(2,2) = 12.0d0 * e%material%modulus * e%moment_of_inertia / (l**3)
        ans(3,3) = 4.0d0 * e%material%modulus * e%moment_of_inertia / l
        ans(2,3) = 6.0d0 * e%material%modulus * e%moment_of_inertia / (l**2)
        ans(3,2) = ans(2,3)
        ans(1,4) = -ans(1,1)
        ans(4,1) = ans(1,4)
        ans(2,5) = -ans(2,2)
        ans(5,2) = ans(2,5)
        ans(2,6) = ans(2,3)
        ans(6,2) = ans(2,6)
        ans(3,5) = -ans(2,3)
        ans(5,3) = ans(3,5)
        ans(3,6) = 2.0d0 * e%material%modulus * e%moment_of_inertia / l
        ans(6,3) = ans(3,6)
        ans(4:6,4:6) = ans(1:3,1:3)
        ans(5,6) = -ans(2,3)
        ans(6,5) = ans(5,6)

        ans = matmul(transpose(T), matmul(ans, T))
    
        ! Test
        k = e%stiffness_matrix()
        if (.not.assert(ans, k, tol * maxval(abs(ans)))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_stiffness_matrix -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_mass_matrix() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Local Variables
        real(real64) :: l, x1, y1, x2, y2, w, h, f
        real(real64), allocatable, dimension(:,:) :: m, ans, T
        type(beam_element_2d) :: e

        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(y1)
        call random_number(y2)
        e%node_1%index = 1
        e%node_1%x = x1
        e%node_1%y = y1
        e%node_1%z = 0.0d0
        e%node_2%index = 2
        e%node_2%x = x2
        e%node_2%y = y2
        e%node_2%z = 0.0d0
        e%node_1%dof = 3
        e%node_2%dof = 3
        l = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        T = e%rotation_matrix()

        ! Define the cross-sectional properties
        call random_number(w)
        call random_number(h)
        e%area = w * h
        e%moment_of_inertia = w * h**3 / 12.0d0

        ! Define the material properties
        e%material%density = 0.101d0 / 3.86d2
        e%material%modulus = 10.0d6
        e%material%poissons_ratio = 0.33d0

        ! Define the solution
        f = e%material%density * e%area * l / 4.2d2
        allocate(ans(6, 6), source = 0.0d0)
        ans(1,1) = 1.4d2 * f
        ans(4,1) = 7.0d1 * f
        ans(2,2) = 1.56d2 * f
        ans(3,2) = 2.2d1 * l * f
        ans(5,2) = 5.4d1 * f
        ans(6,2) = -1.3d1 * l * f
        ans(2,3) = ans(3,2)
        ans(3,3) = 4.0d0 * l**2 * f
        ans(5,3) = 1.3d1 * l * f
        ans(6,3) = -3.0d0 * l**2 * f
        ans(1,4) = ans(4,1)
        ans(4,4) = ans(1,1)
        ans(2,5) = ans(5,2)
        ans(3,5) = ans(5,3)
        ans(5,5) = ans(2,2)
        ans(6,5) = -2.2d1 * l * f
        ans(2,6) = ans(6,2)
        ans(3,6) = ans(6,3)
        ans(5,6) = ans(6,5)
        ans(6,6) = ans(3,3)
        ans = matmul(transpose(T), matmul(ans, T))

        ! Test
        m = e%mass_matrix()
        if (.not.assert(ans, m, tol * maxval(abs(ans)))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_mass_matrix -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_beam2d_ext_force() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6
    
        ! Local Variables
        real(real64) :: l, x1, y1, x2, y2, w, h
        real(real64) :: F(6), ans(6), q(2), Fe(6,2)
        real(real64), allocatable, dimension(:,:) :: T
        type(beam_element_2d) :: e
    
        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(y1)
        call random_number(y2)
        e%node_1%index = 1
        e%node_1%x = x1
        e%node_1%y = y1
        e%node_1%z = 0.0d0
        e%node_2%index = 2
        e%node_2%x = x2
        e%node_2%y = y2
        e%node_2%z = 0.0d0
        e%node_1%dof = 3
        e%node_2%dof = 3
        l = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        T = e%rotation_matrix()
        q = [0.0d0, 1.0d0]
    
        ! Define the cross-sectional properties
        call random_number(w)
        call random_number(h)
        e%area = w * h
        e%moment_of_inertia = w * h**3 / 12.0d0
    
        ! Define the material properties
        e%material%density = 0.101d0 / 3.86d2
        e%material%modulus = 10.0d6
        e%material%poissons_ratio = 0.33d0
    
        ! Define the solution
        Fe = l * reshape( &
            [0.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, &
            0.0d0, 0.5d0, l / 1.2d1, 0.0d0, 0.5d0, -l / 1.2d1], [6, 2])
        ans = matmul(Fe, q)
        ans = matmul(T, ans)

        ! Test
        F = e%external_force_vector(q)
        if (.not.assert(ans, F, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_beam2d_ext_force -1"
        end if
    end function
    
! ------------------------------------------------------------------------------
    function test_boundary_conditions() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 50
        integer(int32), parameter :: bc_index = 25
        real(real64), parameter :: bc_value = 0.1d0
        real(real64) :: k(n,n), f(n), check(n), kcheck(n-1,n-1), fcheck(n-1)
        real(real64), allocatable, dimension(:) :: fnew, frestore
        real(real64), allocatable, dimension(:,:) :: knew
        integer(int32) :: gdofs(1)

        ! Initialization
        rst = .true.
        call random_number(k)
        call random_number(f)
        check = 0.0d0
        check(bc_index) = 1.0d0
        kcheck(1:bc_index-1,1:bc_index-1) = k(1:bc_index-1,1:bc_index-1)
        kcheck(1:bc_index-1,bc_index:) = k(1:bc_index-1,bc_index+1:)
        kcheck(bc_index:,1:bc_index-1) = k(bc_index+1:,1:bc_index-1)
        kcheck(bc_index:,bc_index:) = k(bc_index+1:,bc_index+1:)
        gdofs(1) = bc_index
        fcheck(1:bc_index-1) = f(1:bc_index-1)
        fcheck(bc_index:) = f(bc_index+1:)

        ! Start by applying a known displacement condition
        call apply_displacement_constraint(bc_index, bc_value, k, f)

        ! Test the matrix
        if (.not.assert(k(bc_index,:), check)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions -1"
        end if
        if (.not.assert(f(bc_index), bc_value)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions -2"
        end if

        ! Remove the row and column from the matrix
        knew = apply_boundary_conditions(gdofs, k)

        ! Test
        if (.not.assert(kcheck, knew)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions -3"
        end if

        ! Remove the item from the vector
        fnew = apply_boundary_conditions(gdofs, f)

        ! Test
        if (.not.assert(fcheck, fnew)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions -4"
        end if

        ! Restore F
        frestore = restore_constrained_values(gdofs, fnew)
        f(gdofs) = 0.0d0    ! Just for testing
        
        ! Test
        if (.not.assert(frestore, f)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions -5"
        end if
    end function

    ! ------------------------------------------------------------------------------
    function test_boundary_conditions_2() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 50
        integer(int32), parameter :: bc_index = n
        real(real64), parameter :: bc_value = 0.1d0
        real(real64) :: k(n,n), f(n), check(n), kcheck(n-1,n-1), fcheck(n-1)
        real(real64), allocatable, dimension(:) :: fnew, frestore
        real(real64), allocatable, dimension(:,:) :: knew
        integer(int32) :: gdofs(1)

        ! Initialization
        rst = .true.
        call random_number(k)
        call random_number(f)
        check = 0.0d0
        check(bc_index) = 1.0d0
        kcheck(1:bc_index-1,1:bc_index-1) = k(1:bc_index-1,1:bc_index-1)
        kcheck(1:bc_index-1,bc_index:) = k(1:bc_index-1,bc_index+1:)
        kcheck(bc_index:,1:bc_index-1) = k(bc_index+1:,1:bc_index-1)
        kcheck(bc_index:,bc_index:) = k(bc_index+1:,bc_index+1:)
        gdofs(1) = bc_index
        fcheck(1:bc_index-1) = f(1:bc_index-1)
        fcheck(bc_index:) = f(bc_index+1:)

        ! Start by applying a known displacement condition
        call apply_displacement_constraint(bc_index, bc_value, k, f)

        ! Test the matrix
        if (.not.assert(k(bc_index,:), check)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_2 -1"
        end if
        if (.not.assert(f(bc_index), bc_value)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_2 -2"
        end if

        ! Remove the row and column from the matrix
        knew = apply_boundary_conditions(gdofs, k)

        ! Test
        if (.not.assert(kcheck, knew)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_2 -3"
        end if

        ! Remove the item from the vector
        fnew = apply_boundary_conditions(gdofs, f)

        ! Test
        if (.not.assert(fcheck, fnew)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_2 -4"
        end if

        ! Restore F
        frestore = restore_constrained_values(gdofs, fnew)
        f(gdofs) = 0.0d0    ! Just for testing
        
        ! Test
        if (.not.assert(frestore, f)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_2 -5"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_boundary_conditions_csr() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: n = 20
        integer(int32), parameter :: nbc = 3
        real(real64) :: k(n,n), f(n), fcsr(n)
        real(real64), allocatable, dimension(:,:) :: knew_dense, knew_csr, &
            restored_dense, restored_csr
        type(csr_matrix) :: kcsr, knewcsr, restoredcsr
        integer(int32) :: gdofs_dense(nbc), gdofs_csr(nbc)
        integer(int32) :: free_dofs(n - nbc)
        integer(int32) :: i, j

        ! Initialization
        rst = .true.
        call random_number(k)
        call random_number(f)
        gdofs_dense = [3, 10, 17]
        gdofs_csr = gdofs_dense
        kcsr = k

        ! Remove the rows and columns via the dense and CSR implementations
        knew_dense = apply_boundary_conditions(gdofs_dense, k)
        knewcsr = apply_boundary_conditions(gdofs_csr, kcsr)
        allocate(knew_csr(n - nbc, n - nbc))
        knew_csr = knewcsr

        ! Test
        if (.not.assert(knew_dense, knew_csr)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_csr -1"
        end if

        ! Restore the reduced matrices and compare the dense and CSR paths.
        allocate(restored_dense(n,n), restored_csr(n,n))
        restored_dense = 0.0d0
        j = 0
        do i = 1, n
            if (i /= gdofs_dense(1) .and. i /= gdofs_dense(2) .and. &
                i /= gdofs_dense(3)) then
                j = j + 1
                free_dofs(j) = i
            end if
        end do
        restored_dense(free_dofs, free_dofs) = knew_dense
        restoredcsr = restore_constrained_values(gdofs_csr, knewcsr)
        restored_csr = restoredcsr
        if (.not.assert(restored_dense, restored_csr)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_csr -2"
        end if

        ! Compare dense and CSR displacement constraints.
        fcsr = f
        kcsr = k
        call apply_displacement_constraint(gdofs_dense(2), 0.25d0, k, f)
        call apply_displacement_constraint(gdofs_csr(2), 0.25d0, kcsr, fcsr)
        restored_csr = kcsr
        if (.not.assert(k, restored_csr) .or. &
            .not.assert(f, fcsr)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_boundary_conditions_csr -3"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_connectivity_matrix() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
        integer(int32), parameter :: gdof = 9 ! 3 dof per node, 3 nodes
        real(real64) :: x1, y1, x2, y2, x3, y3, a1, a2, i1, i2, &
            L1ans(6,9), L2ans(6,9), L1(6,9), L2(6,9)
        type(csr_matrix) :: L1csr, L2csr
        type(beam_element_2d) :: b1, b2
        type(material) :: mat
        type(node) :: nodes(3)

        ! Initialization
        rst = .true.
        call random_number(x1)
        call random_number(x2)
        call random_number(x3)
        call random_number(y1)
        call random_number(y2)
        call random_number(y3)
        call random_number(a1)
        call random_number(a2)
        call random_number(i1)
        call random_number(i2)
        mat%modulus = 1.0d7
        mat%density = 0.1d0 / 3.86d2
        mat%poissons_ratio = 0.33d0
        L1ans = 0.0d0
        L2ans = 0.0d0

        L1ans(1,1) = 1.0d0
        L1ans(2,2) = 1.0d0
        L1ans(3,3) = 1.0d0
        L1ans(4,4) = 1.0d0
        L1ans(5,5) = 1.0d0
        L1ans(6,6) = 1.0d0

        L2ans(1,4) = 1.0d0
        L2ans(2,5) = 1.0d0
        L2ans(3,6) = 1.0d0
        L2ans(4,7) = 1.0d0
        L2ans(5,8) = 1.0d0
        L2ans(6,9) = 1.0d0

        ! Ensure x1 < x2 < x3
        x2 = x2 + x1
        x3 = x2 + x3

        ! Create a mesh of two beam elements as follows:
        !
        ! 1 ----- 2 ----- 3
        b1%node_1%index = 1
        b1%node_1%x = x1
        b1%node_1%y = y1
        b1%node_1%z = 0.0d0
        b1%node_1%dof = 3
        b1%node_2%index = 2
        b1%node_2%x = x2
        b1%node_2%y = y2
        b1%node_2%z = 0.0d0
        b1%node_2%dof = 3
        b1%area = a1
        b1%moment_of_inertia = i1
        b1%material = mat

        b2%node_1%index = 2
        b2%node_1%x = x2
        b2%node_1%y = y2
        b2%node_1%z = 0.0d0
        b2%node_1%dof = 3
        b2%node_2%index = 3
        b2%node_2%x = x3
        b2%node_2%y = y3
        b2%node_2%z = 0.0d0
        b2%node_2%dof = 3
        b2%area = a2
        b2%moment_of_inertia = i2
        b2%material = mat

        ! Initialize the node list
        nodes = [b1%node_1, b1%node_2, b2%node_2]

        ! Construct the matrices
        L1csr = create_connectivity_matrix(gdof, b1, nodes)
        L2csr = create_connectivity_matrix(gdof, b2, nodes)
        L1 = L1csr
        L2 = L2csr

        ! Tests
        if (.not.assert(L1, L1ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_connectivity_matrix -1"
        end if
        if (.not.assert(L2, L2ans)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_connectivity_matrix -2"
        end if
    end function

! ------------------------------------------------------------------------------
function test_global_assembly() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer(int32) :: i, j, eidx, offset
    real(real64) :: local_k(6,6), local_m(6,6)
    real(real64) :: expected_k(9,9), expected_m(9,9)
    real(real64), allocatable :: actual_k(:,:), actual_m(:,:)
    type(beam_element_2d) :: elements(2)
    type(material) :: mat
    type(node) :: nodes(3)
    type(csr_matrix) :: kcsr, mcsr

    ! Initialization
    rst = .true.
    mat = material(2.0d0, 10.0d0, 0.25d0)
    nodes = [ &
        node(1, 3, 0.0d0, 0.0d0, 0.0d0), &
        node(2, 3, 1.0d0, 0.0d0, 0.0d0), &
        node(3, 3, 2.0d0, 0.0d0, 0.0d0) &
    ]
    elements(1) = beam_element_2d(mat, 1.0d0, 0.5d0, nodes(1), nodes(2))
    elements(2) = beam_element_2d(mat, 1.0d0, 0.5d0, nodes(2), nodes(3))
    expected_k = 0.0d0
    expected_m = 0.0d0

    do eidx = 1, 2
        if (eidx == 1) then
            offset = 1
        else
            offset = 4
        end if
        local_k = elements(eidx)%stiffness_matrix()
        local_m = elements(eidx)%mass_matrix()
        do i = 1, 6
            do j = 1, 6
                expected_k(offset + i - 1, offset + j - 1) = &
                    expected_k(offset + i - 1, offset + j - 1) + local_k(i,j)
                expected_m(offset + i - 1, offset + j - 1) = &
                    expected_m(offset + i - 1, offset + j - 1) + local_m(i,j)
            end do
        end do
    end do

    call assemble_static_system(9, elements, nodes, kcsr)
    allocate(actual_k(9,9))
    actual_k = kcsr
    if (.not.assert(actual_k, expected_k)) then
        rst = .false.
    end if
    call assemble_static_system(9, elements, nodes, actual_k)
    if (.not.assert(actual_k, expected_k)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_global_assembly -3"
    end if
    if (.not.rst) print "(A)", "TEST FAILED: test_global_assembly -1"

    call assemble_dynamic_system(9, elements, nodes, mcsr, kcsr)
    deallocate(actual_k)
    allocate(actual_m(9,9), actual_k(9,9))
    actual_m = mcsr
    actual_k = kcsr
    if (.not.assert(actual_m, expected_m) .or. &
        .not.assert(actual_k, expected_k)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_global_assembly -2"
    end if
    call assemble_dynamic_system(9, elements, nodes, actual_m, actual_k)
    if (.not.assert(actual_m, expected_m) .or. &
        .not.assert(actual_k, expected_k)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_global_assembly -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_shape_function_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6


    ! Local Variables
    real(real64) :: x1, y1, z1, x2, y2, z2, A, rho, E, nu, Ixx, Iyy, Izz, L, G, s
    real(real64), allocatable, dimension(:,:) :: N, ans
    type(beam_element_3d) :: b
    type(material) :: mat

    ! Initialization
    rst = .true.
    E = 1.0d7
    nu = 0.33d0
    rho = 0.1d0 / 3.86d2
    G = E / (2.0d0 * (1.0d0 + nu))
    call random_number(x1)
    call random_number(y1)
    call random_number(z1)
    call random_number(x2)
    call random_number(y2)
    call random_number(z2)
    call random_number(A)
    call random_number(Ixx)
    call random_number(Iyy)
    call random_number(Izz)
    b%area = A
    b%Ixx = Ixx
    b%Iyy = Iyy
    b%Izz = Izz
    b%material%density = rho
    b%material%modulus = E
    b%material%poissons_ratio = nu
    b%node_1%dof = 6
    b%node_1%index = 1
    b%node_1%x = x1
    b%node_1%y = y1
    b%node_1%z = z1
    b%node_2%dof = 6
    b%node_2%index = 2
    b%node_2%x = x2 + 1.25d0
    b%node_2%y = y2
    b%node_2%z = z2
    L = b%length()

    ! Define the answer
    call random_number(s)
    allocate(ans(4,12), source = 0.0d0)
    ans(1,1) = 0.5d0 * (1.0d0 - s)
    ans(2,2) = 0.25d0 * (s - 1.0d0)**2 * (s + 2.0d0)
    ans(3,3) = 0.25d0 * (s - 1.0d0)**2 * (s + 2.0d0)
    ans(4,4) = 0.5d0 * (1.0d0 - s)
    ans(3,5) = -0.125d0 * L * (s - 1.0d0)**2 * (s + 1.0d0)
    ans(2,6) = 0.125d0 * L * (s - 1.0d0)**2 * (s + 1.0d0)
    ans(1,7) = 0.5d0 * (1.0d0 + s)
    ans(2,8) = 0.25d0 * (2.0d0 - s) * (s + 1.0d0)**2
    ans(3,9) = 0.25d0 * (2.0d0 - s) * (s + 1.0d0)**2
    ans(4,10) = 0.5d0 * (1.0d0 + s)
    ans(3,11) = 0.125d0 * L * (1.0d0 - s) * (1.0d0 + s)**2
    ans(2,12) = 0.125d0 * L * (s - 1.0d0) * (1.0d0 + s)**2

    ! Compute the matrix
    N = b%shape_function_matrix([s])

    ! Test
    if (.not.assert(N, ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_shape_function_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_strain_displacement() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: x1, y1, z1, x2, y2, z2, A, rho, E, nu, Ixx, Iyy, Izz, L, G, s
    real(real64), allocatable, dimension(:,:) :: Be, ans
    type(beam_element_3d) :: b
    type(material) :: mat

    ! Initialization
    rst = .true.
    E = 1.0d7
    nu = 0.33d0
    rho = 0.1d0 / 3.86d2
    G = E / (2.0d0 * (1.0d0 + nu))
    call random_number(x1)
    call random_number(y1)
    call random_number(z1)
    call random_number(x2)
    call random_number(y2)
    call random_number(z2)
    call random_number(A)
    call random_number(Ixx)
    call random_number(Iyy)
    call random_number(Izz)
    b%area = A
    b%Ixx = Ixx
    b%Iyy = Iyy
    b%Izz = Izz
    b%material%density = rho
    b%material%modulus = E
    b%material%poissons_ratio = nu
    b%node_1%dof = 6
    b%node_1%index = 1
    b%node_1%x = x1
    b%node_1%y = y1
    b%node_1%z = z1
    b%node_2%dof = 6
    b%node_2%index = 2
    b%node_2%x = x2 + 1.25d0
    b%node_2%y = y2
    b%node_2%z = z2
    L = b%length()

    ! Define the answer
    call random_number(s)
    allocate(ans(4,12), source = 0.0d0)
    ans(1,1) = -1.0d0 / L
    ans(2,2) = 6.0d0 * s / L**2
    ans(3,3) = 6.0d0 * s / L**2
    ans(4,4) = -1.0d0 / L
    ans(3,5) = (1.0d0 - 3.0d0 * s) / L
    ans(2,6) = (3.0d0 * s - 1.0d0) / L
    ans(1,7) = 1.0d0 / L
    ans(2,8) = -6.0d0 * s / L**2
    ans(3,9) = -6.0d0 * s / L**2
    ans(4,10) = 1.0d0 / L
    ans(3,11) = -(3.0d0 * s + 1.0d0) / L
    ans(2,12) = (3.0d0 * s + 1.0d0) / L

    ! Compute the matrix
    Be = b%strain_displacement_matrix([s])

    ! Test
    if (.not.assert(Be, ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_strain_displacement -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_stiffness_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, A, rho, E, nu, Ixx, &
        Iyy, Izz, Iyz, L, G, s
    real(real64), allocatable, dimension(:,:) :: T, K, ans
    type(beam_element_3d) :: b
    type(material) :: mat

    ! Initialization
    rst = .true.
    E = 1.0d7
    nu = 0.33d0
    rho = 0.1d0 / 3.86d2
    G = E / (2.0d0 * (1.0d0 + nu))
    call random_number(x1)
    call random_number(y1)
    call random_number(z1)
    call random_number(x2)
    call random_number(y2)
    call random_number(z2)
    call random_number(x3)
    call random_number(y3)
    call random_number(z3)
    call random_number(A)
    call random_number(Ixx)
    call random_number(Iyy)
    call random_number(Izz)
    call random_number(Iyz)
    b%area = A
    b%Ixx = Ixx
    b%Iyy = Iyy
    b%Izz = Izz
    b%Iyz = Iyz
    b%material%density = rho
    b%material%modulus = E
    b%material%poissons_ratio = nu
    b%node_1%dof = 6
    b%node_1%index = 1
    b%node_1%x = x1
    b%node_1%y = y1
    b%node_1%z = z1
    b%node_2%dof = 6
    b%node_2%index = 2
    b%node_2%x = x2
    b%node_2%y = y2
    b%node_2%z = z2
    b%orientation_point%x = x3
    b%orientation_point%y = y3
    b%orientation_point%z = z3
    L = b%length()
    T = b%rotation_matrix()

    ! Define the answer
    allocate(ans(12, 12), source = 0.0d0)
    ans(1,1) = A * E / L
    ans(7,1) = -ans(1,1)
    ans(2,2) = 1.2d1 * E * Izz / L**3
    ans(6,2) = 6.0d0 * E * Izz / L**2
    ans(8,2) = -ans(2,2)
    ans(12,2) = ans(6,2)
    ans(3,3) = 1.2d1 * E * Iyy / L**3
    ans(5,3) = -6.0d0 * E * Iyy / L**2
    ans(9,3) = -ans(3,3)
    ans(11,3) = ans(5,3)
    ans(2,3) = 1.2d1 * E * Iyz / L**3
    ans(5,2) = -6.0d0 * E * Iyz / L**2
    ans(8,3) = -ans(2,3)
    ans(11,2) = 6.0d0 * E * Iyz / L**2
    ans(6,3) = 6.0d0 * E * Iyz / L**2
    ans(8,5) = 6.0d0 * E * Iyz / L**2
    ans(12,3) = -6.0d0 * E * Iyz / L**2
    ans(5,9) = 6.0d0 * E * Iyz / L**2
    ans(11,9) = -6.0d0 * E * Iyz / L**2
    ans(6,5) = -4.0d0 * E * Iyz / L
    ans(12,5) = -2.0d0 * E * Iyz / L
    ans(6,11) = -2.0d0 * E * Iyz / L
    ans(12,11) = -4.0d0 * E * Iyz / L
    ans(4,4) = G * Ixx / L
    ans(10,4) = -ans(4,4)
    ans(3,5) = ans(5,3)
    ans(2,5) = ans(5,2)
    ans(2,9) = -ans(2,3)
    ans(2,11) = ans(11,2)
    ans(3,6) = ans(6,3)
    ans(3,8) = ans(8,3)
    ans(3,12) = ans(12,3)
    ans(5,5) = 4.0d0 * E * Iyy / L
    ans(9,5) = -ans(11,3)
    ans(11,5) = 2.0d0 * E * Iyy / L
    ans(2,6) = ans(6,2)
    ans(6,6) = 4.0d0 * E * Izz / L
    ans(8,6) = -ans(12,2)
    ans(12,6) = 2.0d0 * E * Izz / L
    ans(1,7) = ans(7,1)
    ans(7,7) = ans(1,1)
    ans(2,8) = ans(8,2)
    ans(5,6) = ans(6,5)
    ans(5,8) = ans(8,5)
    ans(5,12) = ans(12,5)
    ans(6,8) = ans(8,6)
    ans(8,8) = ans(2,2)
    ans(12,8) = -ans(6,2)
    ans(3,9) = ans(9,3)
    ans(6,9) = -6.0d0 * E * Iyz / L**2
    ans(8,9) = -ans(2,9)
    ans(9,2) = ans(2,9)
    ans(9,6) = ans(6,9)
    ans(9,8) = ans(8,9)
    ans(5,9) = ans(9,5)
    ans(9,9) = ans(3,3)
    ans(11,9) = -ans(5,3)
    ans(4,10) = ans(10,4)
    ans(10,10) = ans(4,4)
    ans(3,11) = ans(11,3)
    ans(8,11) = -6.0d0 * E * Iyz / L**2
    ans(9,11) = ans(11,9)
    ans(11,6) = ans(6,11)
    ans(11,8) = ans(8,11)
    ans(5,11) = ans(11,5)
    ans(9,11) = ans(11,9)
    ans(11,11) = ans(5,5)
    ans(2,12) = ans(12,2)
    ans(9,12) = -6.0d0 * E * Iyz / L**2
    ans(12,9) = ans(9,12)
    ans(11,12) = ans(12,11)
    ans(6,12) = ans(12,6)
    ans(8,12) = ans(12,8)
    ans(12,12) = ans(6,6)
    ans = matmul(transpose(T), matmul(ans, T))

    ! Compute the matrix
    K = b%stiffness_matrix()

    ! Test
    if (.not.assert(K, ans, tol * maxval(abs(ans)))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_stiffness_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_bending_stress() result(rst)
    logical :: rst

    real(real64), parameter :: length = 5.0d0
    real(real64), parameter :: area = 2.5d0
    real(real64), parameter :: ixx = 0.8d0
    real(real64), parameter :: iyy = 1.2d0
    real(real64), parameter :: izz = 1.6d0
    real(real64), parameter :: iyz = -0.3d0
    real(real64), parameter :: modulus = 2.0d7
    real(real64), parameter :: ky = 2.0d-3
    real(real64), parameter :: kz = -1.0d-3
    real(real64), parameter :: tol = 1.0d-12
    real(real64), parameter :: s(1) = [0.25d0]

    real(real64) :: displacement(12), expected(4)
    real(real64), allocatable :: stress(:)
    type(material) :: mat
    type(beam_element_3d) :: e

    rst = .true.
    mat = material(modulus, 0.25d0, 1.0d0)
    e = beam_element_3d(mat, area, ixx, iyy, izz, iyz, &
        node(1, 6, 0.0d0, 0.0d0, 0.0d0), &
        node(2, 6, length, 0.0d0, 0.0d0), point(0.0d0, 0.0d0, -1.0d0))
    displacement = 0.0d0
    displacement(8) = 0.5d0 * ky * length**2
    displacement(12) = ky * length
    displacement(9) = 0.5d0 * kz * length**2
    displacement(11) = -kz * length
    expected = [0.0d0, modulus * (izz * ky + iyz * kz), &
        modulus * (iyz * ky + iyy * kz), 0.0d0]
    stress = e%stress(displacement, s)

    if (.not.assert(stress, expected, tol * maxval(abs(expected)))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_bending_stress -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_bending_strain() result(rst)
    logical :: rst

    real(real64), parameter :: length = 5.0d0
    real(real64), parameter :: ky = 2.0d-3
    real(real64), parameter :: kz = -1.0d-3
    real(real64), parameter :: s(1) = [0.25d0]

    real(real64) :: displacement(12), expected(4)
    real(real64), allocatable :: strain(:)
    type(material) :: mat
    type(beam_element_3d) :: e

    rst = .true.
    mat = material(2.0d7, 0.25d0, 1.0d0)
    e = beam_element_3d(mat, 2.5d0, 0.8d0, 1.2d0, 1.6d0, -0.3d0, &
        node(1, 6, 0.0d0, 0.0d0, 0.0d0), &
        node(2, 6, length, 0.0d0, 0.0d0), point(0.0d0, 0.0d0, -1.0d0))
    displacement = 0.0d0
    displacement(8) = 0.5d0 * ky * length**2
    displacement(12) = ky * length
    displacement(9) = 0.5d0 * kz * length**2
    displacement(11) = -kz * length
    expected = [0.0d0, ky, kz, 0.0d0]
    strain = e%strain(displacement, s)

    if (.not.assert(strain, expected, 1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_bending_strain -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_internal_results() result(rst)
    logical :: rst

    real(real64), parameter :: length = 5.0d0
    real(real64), parameter :: ixx = 0.8d0
    real(real64), parameter :: iyy = 1.2d0
    real(real64), parameter :: izz = 1.6d0
    real(real64), parameter :: iyz = -0.3d0
    real(real64), parameter :: modulus = 2.0d7
    real(real64), parameter :: cy = 1.0d-4
    real(real64), parameter :: cz = -2.0d-4
    real(real64), parameter :: s(1) = [0.25d0]
    real(real64), parameter :: x = 0.5d0 * length * (s(1) + 1.0d0)
    real(real64), parameter :: tol = 1.0d-12

    real(real64) :: displacement(12), expected_moment(3), expected_shear(2)
    type(material) :: mat
    type(beam_element_3d) :: e

    rst = .true.
    mat = material(modulus, 0.25d0, 1.0d0)
    e = beam_element_3d(mat, 2.5d0, ixx, iyy, izz, iyz, &
        node(1, 6, 0.0d0, 0.0d0, 0.0d0), &
        node(2, 6, length, 0.0d0, 0.0d0), point(0.0d0, 0.0d0, -1.0d0))
    displacement = 0.0d0
    displacement(8) = cy * length**3
    displacement(12) = 3.0d0 * cy * length**2
    displacement(9) = cz * length**3
    displacement(11) = -3.0d0 * cz * length**2
    expected_moment = [0.0d0, 0.0d0, &
        modulus * (izz * 6.0d0 * cy * x + iyz * 6.0d0 * cz * x)]
    expected_moment(2) = modulus * (iyz * 6.0d0 * cy * x + &
        iyy * 6.0d0 * cz * x)
    expected_shear = [modulus * (izz * 6.0d0 * cy + iyz * 6.0d0 * cz), &
        -modulus * (iyz * 6.0d0 * cy + iyy * 6.0d0 * cz)]

    if (.not.assert(e%bending_moment(displacement, s), expected_moment, &
            tol * maxval(abs(expected_moment)))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_internal_results -1"
    end if
    if (.not.assert(e%shear_force(displacement, s), expected_shear, &
            tol * maxval(abs(expected_shear)))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_internal_results -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_constitutive_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    real(real64) :: area, E, nu, Ixx, Iyy, Izz, Iyz, G
    real(real64), allocatable, dimension(:,:) :: D, ans
    type(beam_element_3d) :: b

    ! Initialization
    rst = .true.
    area = 2.5d0
    E = 2.0d7
    nu = 0.25d0
    Ixx = 0.8d0
    Iyy = 1.2d0
    Izz = 1.6d0
    Iyz = -0.3d0
    G = E / (2.0d0 * (1.0d0 + nu))
    b%area = area
    b%Ixx = Ixx
    b%Iyy = Iyy
    b%Izz = Izz
    b%Iyz = Iyz
    b%material%modulus = E
    b%material%poissons_ratio = nu

    ! Define the expected constitutive matrix, including bending coupling.
    allocate(ans(4,4), source = 0.0d0)
    ans(1,1) = area * E
    ans(2,2) = Izz * E
    ans(2,3) = Iyz * E
    ans(3,2) = ans(2,3)
    ans(3,3) = Iyy * E
    ans(4,4) = Ixx * G

    ! Test
    D = b%constitutive_matrix()
    if (.not.assert(D, ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_constitutive_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_beam3d_mass_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, A, rho, E, nu, Ixx, &
        Iyy, Izz, L, G, s, v(12)
    real(real64), allocatable, dimension(:,:) :: T, M, ans
    type(beam_element_3d) :: b
    type(material) :: mat

    ! Initialization
    rst = .true.
    E = 1.0d7
    nu = 0.33d0
    rho = 0.1d0 / 3.86d2
    G = E / (2.0d0 * (1.0d0 + nu))
    call random_number(x1)
    call random_number(y1)
    call random_number(z1)
    call random_number(x2)
    call random_number(y2)
    call random_number(z2)
    call random_number(x3)
    call random_number(y3)
    call random_number(z3)
    call random_number(A)
    call random_number(Ixx)
    call random_number(Iyy)
    call random_number(Izz)
    b%area = A
    b%Ixx = Ixx
    b%Iyy = Iyy
    b%Izz = Izz
    b%material%density = rho
    b%material%modulus = E
    b%material%poissons_ratio = nu
    b%node_1%dof = 6
    b%node_1%index = 1
    b%node_1%x = x1
    b%node_1%y = y1
    b%node_1%z = z1
    b%node_2%dof = 6
    b%node_2%index = 2
    b%node_2%x = x2
    b%node_2%y = y2
    b%node_2%z = z2
    b%orientation_point%x = x3
    b%orientation_point%y = y3
    b%orientation_point%z = z3
    L = b%length()
    T = b%rotation_matrix()

    ! Define the answer (m = rho A L; torsion uses rho Ixx L)
    s = rho * A * L
    allocate(ans(12, 12), source = 0.0d0)
    ans(1,1) = s / 3.0d0
    ans(7,1) = s / 6.0d0
    ans(2,2) = 1.3d1 * s / 3.5d1
    ans(6,2) = 1.1d1 * s * L / 2.1d2
    ans(8,2) = 9.0d0 * s / 7.0d1
    ans(12,2) = -1.3d1 * s * L / 4.2d2
    ans(3,3) = 1.3d1 * s / 3.5d1
    ans(5,3) = -1.1d1 * s * L / 2.1d2
    ans(9,3) = 9.0d0 * s / 7.0d1
    ans(11,3) = 1.3d1 * s * L / 4.2d2
    ans(4,4) = rho * Ixx * L / 3.0d0
    ans(10,4) = rho * Ixx * L / 6.0d0
    ans(3,5) = ans(5,3)
    ans(5,5) = s * L**2 / 1.05d2
    ans(9,5) = -1.3d1 * s * L / 4.2d2
    ans(11,5) = -s * L**2 / 1.4d2
    ans(2,6) = ans(6,2)
    ans(6,6) = s * L**2 / 1.05d2
    ans(8,6) = 1.3d1 * s * L / 4.2d2
    ans(12,6) = -s * L**2 / 1.4d2
    ans(1,7) = ans(7,1)
    ans(7,7) = ans(1,1)
    ans(2,8) = ans(8,2)
    ans(6,8) = ans(8,6)
    ans(8,8) = ans(2,2)
    ans(12,8) = -1.1d1 * s * L / 2.1d2
    ans(3,9) = ans(9,3)
    ans(5,9) = ans(9,5)
    ans(9,9) = ans(3,3)
    ans(11,9) = 1.1d1 * s * L / 2.1d2
    ans(4,10) = ans(10,4)
    ans(10,10) = ans(4,4)
    ans(3,11) = ans(11,3)
    ans(5,11) = ans(11,5)
    ans(9,11) = ans(11,9)
    ans(11,11) = ans(5,5)
    ans(2,12) = ans(12,2)
    ans(6,12) = ans(12,6)
    ans(8,12) = ans(12,8)
    ans(12,12) = ans(6,6)
    ans = matmul(transpose(T), matmul(ans, T))

    ! Compute the matrix
    M = b%mass_matrix()

    ! Test
    if (.not.assert(M, ans, tol * maxval(abs(ans)))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_mass_matrix -1"
    end if

    ! A rigid translation must carry the total element mass
    v = 0.0d0
    v([1, 7]) = 1.0d0
    if (.not.is_symmetric(M) .or. &
        abs(dot_product(v, matmul(M, v)) - rho * A * L) > &
        tol * rho * A * L) then
        rst = .false.
        print "(A)", "TEST FAILED: test_beam3d_mass_matrix -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_truss_elements() result(rst)
    logical :: rst
    integer(int32) :: row, col
    real(real64), parameter :: tol = 1.0d-12
    real(real64), dimension(4) :: axial_2d
    real(real64), dimension(6) :: axial_3d
    real(real64), dimension(4,4) :: expected_k2, expected_m2
    real(real64), dimension(6,6) :: expected_k3, expected_m3
    real(real64), allocatable, dimension(:) :: strain_result
    real(real64), allocatable, dimension(:,:) :: mass, stiffness
    type(material) :: mat
    type(node), dimension(2) :: nodes_2d, nodes_3d
    type(truss_element_2d), dimension(1) :: bars_2d
    type(truss_element_3d), dimension(1) :: bars_3d

    rst = .true.
    mat = material(100.0d0, 0.3d0, 6.0d0)
    nodes_2d(1) = node(1, 2, 0.0d0, 0.0d0, 0.0d0)
    nodes_2d(2) = node(2, 2, 3.0d0, 4.0d0, 0.0d0)
    bars_2d(1) = truss_element_2d(mat, 0.02d0, nodes_2d(1), nodes_2d(2))
    axial_2d = [-0.6d0, -0.8d0, 0.6d0, 0.8d0]
    expected_k2 = 0.0d0
    expected_m2 = 0.0d0
    do col = 1, 4
        do row = 1, 4
            expected_k2(row,col) = 0.4d0 * axial_2d(row) * axial_2d(col)
        end do
        expected_m2(col,col) = 0.2d0
    end do
    expected_m2(1,3) = 0.1d0
    expected_m2(3,1) = 0.1d0
    expected_m2(2,4) = 0.1d0
    expected_m2(4,2) = 0.1d0
    call assemble_dynamic_system(4, bars_2d, nodes_2d, mass, stiffness)
    strain_result = bars_2d(1)%strain([0.0d0, 0.0d0, 0.006d0, 0.008d0], [0.0d0])
    if (maxval(abs(stiffness - expected_k2)) > tol .or. &
        maxval(abs(mass - expected_m2)) > tol .or. &
        abs(bars_2d(1)%length() - 5.0d0) > tol .or. &
        abs(strain_result(1) - 0.002d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_truss_elements - 2D"
    end if

    nodes_3d(1) = node(1, 3, 0.0d0, 0.0d0, 0.0d0)
    nodes_3d(2) = node(2, 3, 0.0d0, 0.0d0, 2.0d0)
    bars_3d(1) = truss_element_3d(mat, 0.02d0, nodes_3d(1), nodes_3d(2))
    axial_3d = [0.0d0, 0.0d0, -1.0d0, 0.0d0, 0.0d0, 1.0d0]
    expected_k3 = 0.0d0
    expected_m3 = 0.0d0
    do col = 1, 6
        do row = 1, 6
            expected_k3(row,col) = 1.0d0 * axial_3d(row) * axial_3d(col)
        end do
        expected_m3(col,col) = 0.08d0
    end do
    expected_m3(1,4) = 0.04d0
    expected_m3(4,1) = 0.04d0
    expected_m3(2,5) = 0.04d0
    expected_m3(5,2) = 0.04d0
    expected_m3(3,6) = 0.04d0
    expected_m3(6,3) = 0.04d0
    call assemble_dynamic_system(6, bars_3d, nodes_3d, mass, stiffness)
    if (maxval(abs(stiffness - expected_k3)) > tol .or. &
        maxval(abs(mass - expected_m3)) > tol .or. &
        abs(bars_3d(1)%length() - 2.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_truss_elements - 3D"
    end if

    nodes_3d(2) = node(2, 3, 3.0d0, 4.0d0, 12.0d0)
    bars_3d(1) = truss_element_3d(mat, 0.02d0, nodes_3d(1), nodes_3d(2))
    axial_3d = [-3.0d0, -4.0d0, -12.0d0, 3.0d0, 4.0d0, 12.0d0] / 13.0d0
    do col = 1, 6
        do row = 1, 6
            expected_k3(row,col) = 2.0d0 / 13.0d0 * axial_3d(row) * axial_3d(col)
        end do
    end do
    strain_result = bars_3d(1)%strain([0.0d0, 0.0d0, 0.0d0, &
        0.003d0, 0.004d0, 0.012d0], [0.0d0])
    call assemble_dynamic_system(6, bars_3d, nodes_3d, mass, stiffness)
    if (maxval(abs(stiffness - expected_k3)) > tol .or. &
        abs(strain_result(1) - 0.001d0) > tol .or. &
        abs(bars_3d(1)%length() - 13.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_truss_elements - oblique 3D"
    end if
end function

! ------------------------------------------------------------------------------
function test_generalized_alpha_integrator() result(rst)
    use linalg, only : dense_to_csr
    logical :: rst
    integer(int32) :: index, step, row
    real(real64) :: dt, rho, alpha_m, alpha_f, gamma, beta
    real(real64), dimension(2) :: displacement, velocity, acceleration
    real(real64), dimension(2) :: sparse_displacement, sparse_velocity, sparse_acceleration
    real(real64), dimension(2) :: previous_displacement, previous_velocity, previous_acceleration
    real(real64), dimension(2) :: force_current, force_next, residual
    real(real64), dimension(2) :: default_displacement, default_velocity, default_acceleration
    real(real64), dimension(2) :: object_displacement, object_velocity, object_acceleration
    real(real64), dimension(2) :: sparse_object_displacement, sparse_object_velocity, sparse_object_acceleration
    real(real64), dimension(2) :: expected_displacement, expected_velocity, expected_acceleration
    real(real64), dimension(2,11) :: force_history
    real(real64), dimension(2,2) :: mass, damping, stiffness
    real(real64), dimension(40,40) :: large_mass, large_damping, large_stiffness
    real(real64), dimension(40) :: large_displacement, large_velocity, large_acceleration
    real(real64), dimension(40) :: large_sparse_displacement, large_sparse_velocity, large_sparse_acceleration
    real(real64), dimension(40) :: large_object_displacement, large_object_velocity, large_object_acceleration
    real(real64), dimension(40) :: large_force_current, large_force_next
    type(csr_matrix) :: sparse_mass, sparse_damping, sparse_stiffness
    type(csr_matrix) :: large_sparse_mass, large_sparse_damping, large_sparse_stiffness
    type(dense_generalized_alpha_integrator) :: dense_integrator, reference_integrator
    type(dense_generalized_alpha_integrator) :: default_integrator, large_dense_integrator
    type(sparse_generalized_alpha_integrator) :: sparse_integrator, sparse_reference_integrator
    type(sparse_generalized_alpha_integrator) :: large_integrator, large_reference_integrator
    class(structural_integrator), allocatable :: polymorphic_integrator

    rst = .true.
    dt = 0.05d0
    mass = reshape([2.0d0, 0.2d0, 0.2d0, 1.5d0], [2, 2])
    damping = reshape([0.4d0, -0.1d0, -0.1d0, 0.3d0], [2, 2])
    stiffness = reshape([8.0d0, -2.0d0, -2.0d0, 5.0d0], [2, 2])
    sparse_mass = dense_to_csr(mass)
    sparse_damping = dense_to_csr(damping)
    sparse_stiffness = dense_to_csr(stiffness)

    do index = 1, 3
        rho = 0.5d0 * real(index - 1, real64)
        alpha_m = (2.0d0 * rho - 1.0d0) / (rho + 1.0d0)
        alpha_f = rho / (rho + 1.0d0)
        gamma = 0.5d0 + alpha_f - alpha_m
        beta = 0.25d0 * (1.0d0 + alpha_f - alpha_m)**2
        displacement = [0.3d0, -0.2d0]
        velocity = [0.1d0, 0.4d0]
        force_current = [1.0d0, -0.5d0]
        acceleration = solve_static_system(mass, force_current - &
            matmul(damping, velocity) - matmul(stiffness, displacement))
        sparse_displacement = displacement
        sparse_velocity = velocity
        sparse_acceleration = acceleration
        object_displacement = displacement
        object_velocity = velocity
        object_acceleration = acceleration
        sparse_object_displacement = displacement
        sparse_object_velocity = velocity
        sparse_object_acceleration = acceleration
        call dense_integrator%initialize(mass, damping, stiffness, rho)
        call sparse_integrator%initialize(sparse_mass, sparse_damping, sparse_stiffness, rho)
        call reference_integrator%initialize(mass, damping, stiffness, rho)
        call sparse_reference_integrator%initialize(sparse_mass, sparse_damping, sparse_stiffness, rho)

        do step = 1, 10
            previous_displacement = displacement
            previous_velocity = velocity
            previous_acceleration = acceleration
            force_next = [1.0d0 + 0.2d0 * step, -0.5d0 + 0.1d0 * step]
            if (index == 3 .and. step == 1) then
                default_displacement = displacement
                default_velocity = velocity
                default_acceleration = acceleration
                call default_integrator%initialize(mass, damping, stiffness)
                call default_integrator%step(force_current, force_next, dt, &
                    default_displacement, default_velocity, default_acceleration)
            end if
            call dense_integrator%step(force_current, force_next, dt, &
                displacement, velocity, acceleration)
            call sparse_integrator%step(force_current, force_next, dt, &
                sparse_displacement, sparse_velocity, sparse_acceleration)
            call reference_integrator%step(force_current, force_next, dt, &
                object_displacement, object_velocity, object_acceleration)
            call sparse_reference_integrator%step(force_current, force_next, dt, &
                sparse_object_displacement, sparse_object_velocity, sparse_object_acceleration)

            residual = matmul(mass, (1.0d0 - alpha_m) * acceleration + &
                alpha_m * previous_acceleration) + &
                matmul(damping, (1.0d0 - alpha_f) * velocity + alpha_f * previous_velocity) + &
                matmul(stiffness, (1.0d0 - alpha_f) * displacement + alpha_f * previous_displacement) - &
                ((1.0d0 - alpha_f) * force_next + alpha_f * force_current)
            if (norm2(residual) > 1.0d-10 .or. &
                norm2(displacement - previous_displacement - dt * previous_velocity - &
                    dt**2 * ((0.5d0 - beta) * previous_acceleration + beta * acceleration)) > 1.0d-10 .or. &
                norm2(velocity - previous_velocity - dt * ((1.0d0 - gamma) * previous_acceleration + &
                    gamma * acceleration)) > 1.0d-10 .or. &
                norm2(displacement - sparse_displacement) > 1.0d-8 .or. &
                norm2(velocity - sparse_velocity) > 1.0d-8 .or. &
                norm2(acceleration - sparse_acceleration) > 1.0d-8 .or. &
                norm2(displacement - object_displacement) > 1.0d-10 .or. &
                norm2(velocity - object_velocity) > 1.0d-10 .or. &
                norm2(acceleration - object_acceleration) > 1.0d-10 .or. &
                norm2(displacement - sparse_object_displacement) > 1.0d-8 .or. &
                norm2(velocity - sparse_object_velocity) > 1.0d-8 .or. &
                norm2(acceleration - sparse_object_acceleration) > 1.0d-8) then
                rst = .false.
                print "(A,I0,A,I0)", "TEST FAILED: test_generalized_alpha_integrator -", index, " step ", step
                return
            end if
            if (index == 3 .and. step == 1) then
                if (norm2(displacement - default_displacement) > 1.0d-12 .or. &
                    norm2(velocity - default_velocity) > 1.0d-12 .or. &
                    norm2(acceleration - default_acceleration) > 1.0d-12) then
                    rst = .false.
                    print "(A)", "TEST FAILED: test_generalized_alpha_integrator default"
                    return
                end if
            end if
            force_current = force_next
        end do
        if (index == 2) then
            expected_displacement = displacement
            expected_velocity = velocity
            expected_acceleration = acceleration
        end if
    end do

    force_history(:,1) = [1.0d0, -0.5d0]
    do step = 1, 10
        force_history(:,step+1) = [1.0d0 + 0.2d0 * step, -0.5d0 + 0.1d0 * step]
    end do
    allocate(dense_generalized_alpha_integrator :: polymorphic_integrator)
    select type (polymorphic_integrator)
    type is (dense_generalized_alpha_integrator)
        call polymorphic_integrator%initialize(mass, damping, stiffness, 0.5d0)
    end select
    object_displacement = [0.3d0, -0.2d0]
    object_velocity = [0.1d0, 0.4d0]
    object_acceleration = solve_static_system(mass, force_history(:,1) - &
        matmul(damping, object_velocity) - matmul(stiffness, object_displacement))
    call polymorphic_integrator%solve(force_history, dt, object_displacement, &
        object_velocity, object_acceleration)
    if (norm2(object_displacement - expected_displacement) > 1.0d-10 .or. &
        norm2(object_velocity - expected_velocity) > 1.0d-10 .or. &
        norm2(object_acceleration - expected_acceleration) > 1.0d-10) then
        rst = .false.
        print "(A)", "TEST FAILED: test_generalized_alpha_integrator dense solve"
    end if

    call sparse_integrator%initialize(sparse_mass, sparse_damping, sparse_stiffness, 0.5d0)
    sparse_object_displacement = [0.3d0, -0.2d0]
    sparse_object_velocity = [0.1d0, 0.4d0]
    sparse_object_acceleration = solve_static_system(mass, force_history(:,1) - &
        matmul(damping, sparse_object_velocity) - matmul(stiffness, sparse_object_displacement))
    call sparse_integrator%solve(force_history, dt, sparse_object_displacement, &
        sparse_object_velocity, sparse_object_acceleration)
    if (norm2(sparse_object_displacement - expected_displacement) > 1.0d-8 .or. &
        norm2(sparse_object_velocity - expected_velocity) > 1.0d-8 .or. &
        norm2(sparse_object_acceleration - expected_acceleration) > 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_generalized_alpha_integrator sparse solve"
    end if

    force_next = [2.0d0, 0.5d0]
    call reference_integrator%initialize(mass, damping, stiffness, 0.5d0)
    call reference_integrator%step(force_current, force_next, 0.07d0, &
        expected_displacement, expected_velocity, expected_acceleration)
    call polymorphic_integrator%step(force_current, force_next, 0.07d0, &
        object_displacement, object_velocity, object_acceleration)
    call sparse_integrator%step(force_current, force_next, 0.07d0, &
        sparse_object_displacement, sparse_object_velocity, sparse_object_acceleration)
    if (norm2(object_displacement - expected_displacement) > 1.0d-10 .or. &
        norm2(sparse_object_displacement - expected_displacement) > 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_generalized_alpha_integrator changed dt"
    end if

    large_mass = 0.0d0
    large_damping = 0.0d0
    large_stiffness = 0.0d0
    do row = 1, 40
        large_mass(row,row) = 2.0d0
        large_stiffness(row,row) = 4.0d0
        if (row < 40) then
            large_stiffness(row,row+1) = -1.0d0
            large_stiffness(row+1,row) = -1.0d0
        end if
    end do
    large_sparse_mass = dense_to_csr(large_mass)
    large_sparse_damping = dense_to_csr(large_damping)
    large_sparse_stiffness = dense_to_csr(large_stiffness)
    large_displacement = 0.0d0
    large_velocity = 0.0d0
    large_acceleration = 0.0d0
    large_sparse_displacement = large_displacement
    large_sparse_velocity = large_velocity
    large_sparse_acceleration = large_acceleration
    large_object_displacement = large_displacement
    large_object_velocity = large_velocity
    large_object_acceleration = large_acceleration
    large_force_current = 0.0d0
    large_force_next = [(0.01d0 * real(row, real64), row = 1, 40)]
    call large_dense_integrator%initialize(large_mass, large_damping, large_stiffness, 0.5d0)
    call large_integrator%initialize(large_sparse_mass, large_sparse_damping, &
        large_sparse_stiffness, 0.5d0)
    call large_reference_integrator%initialize(large_sparse_mass, large_sparse_damping, &
        large_sparse_stiffness, 0.5d0)
    call large_dense_integrator%step(large_force_current, large_force_next, dt, &
        large_displacement, large_velocity, large_acceleration)
    call large_integrator%step(large_force_current, large_force_next, dt, &
        large_sparse_displacement, large_sparse_velocity, large_sparse_acceleration)
    call large_reference_integrator%step(large_force_current, large_force_next, dt, &
        large_object_displacement, large_object_velocity, large_object_acceleration)
    if (norm2(large_displacement - large_sparse_displacement) > 1.0d-8 .or. &
        norm2(large_velocity - large_sparse_velocity) > 1.0d-8 .or. &
        norm2(large_acceleration - large_sparse_acceleration) > 1.0d-8 .or. &
        norm2(large_displacement - large_object_displacement) > 1.0d-8 .or. &
        norm2(large_velocity - large_object_velocity) > 1.0d-8 .or. &
        norm2(large_acceleration - large_object_acceleration) > 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_generalized_alpha_integrator sparse GMRES"
    end if

    large_force_current = large_force_next
    large_force_next = 0.5d0 * large_force_current
    call large_dense_integrator%step(large_force_current, large_force_next, 0.07d0, &
        large_displacement, large_velocity, large_acceleration)
    call large_integrator%step(large_force_current, large_force_next, 0.07d0, &
        large_object_displacement, large_object_velocity, large_object_acceleration)
    if (norm2(large_displacement - large_object_displacement) > 1.0d-8 .or. &
        norm2(large_velocity - large_object_velocity) > 1.0d-8 .or. &
        norm2(large_acceleration - large_object_acceleration) > 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_generalized_alpha_integrator sparse changed dt"
    end if
end function

! ------------------------------------------------------------------------------
function test_integration_rules() result(rst)
    ! A linear truss mass matrix integrates a quadratic, so the 2-, 3-, and
    ! 4-point rules are exact and the 1-point rule samples the midpoint.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-12
    real(real64), parameter :: len = 2.0d0
    real(real64), parameter :: area = 0.5d0
    real(real64), parameter :: rho = 3.0d0
    integer(int32) :: rule
    real(real64) :: mtot
    real(real64), allocatable, dimension(:,:) :: m
    type(truss_element_2d) :: truss

    rst = .true.
    mtot = rho * area * len
    truss = truss_element_2d(material(1.0d0, 0.3d0, rho), area, &
        node(1, 2, 0.0d0, 0.0d0, 0.0d0), node(2, 2, len, 0.0d0, 0.0d0))

    m = truss%mass_matrix(DYN_ONE_POINT_INTEGRATION_RULE)
    if (abs(m(1,1) - 0.25d0 * mtot) > tol .or. &
        abs(m(1,3) - 0.25d0 * mtot) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_integration_rules - rule 1"
    end if

    do rule = DYN_TWO_POINT_INTEGRATION_RULE, DYN_FOUR_POINT_INTEGRATION_RULE
        m = truss%mass_matrix(rule)
        if (abs(m(1,1) - mtot / 3.0d0) > tol .or. &
            abs(m(1,3) - mtot / 6.0d0) > tol) then
            rst = .false.
            print "(A, I0)", "TEST FAILED: test_integration_rules - rule ", rule
        end if
    end do
end function

end module