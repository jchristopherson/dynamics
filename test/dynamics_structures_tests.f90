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
        real(real64), parameter :: tol = 5.0d-6
    
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
        if (.not.assert(ans, m, tol)) then
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
    function test_connectivity_matrix() result(rst)
        ! Arguments
        logical :: rst

        ! Variables
    integer(int32), parameter :: gdof = 9 ! 3 dof per node, 3 nodes
        real(real64) :: x1, y1, x2, y2, x3, y3, a1, a2, i1, i2, &
            L1ans(6,9), L2ans(6,9)
        real(real64), allocatable, dimension(:,:) :: L1, L2
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
        L1 = create_connectivity_matrix(gdof, b1, nodes)
        L2 = create_connectivity_matrix(gdof, b2, nodes)

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
end module