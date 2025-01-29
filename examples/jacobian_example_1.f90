module helper
    use iso_fortran_env
    implicit none

contains

pure function finite_difference(t, x) result(rst)
    ! Computes an estimate of dx/dt via finite differences.
    real(real64), intent(in), dimension(:) :: t, x
    real(real64), allocatable, dimension(:) :: rst

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(t)
    allocate(rst(n))
    rst(1) = (x(2) - x(1)) / (t(2) - t(1))
    do i = 2, n - 1
        rst(i) = (x(i+1) - x(i-1)) / (t(i+1) - t(i-1))
    end do
    rst(n) = (x(n) - x(n-1)) / (t(n) - t(n-1))
end function

end module

program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use helper
    implicit none

    ! The purpose of this example is to illustrate the use of a Jacobian matrix
    ! and compare it to a numerical differentiation.  As such, this example
    ! defines the joint variables and solves the forward problem to determine
    ! the end-effector motion as a function of time.  
    !
    ! This example follows Example 214 out of Jazar's text: "Theory of Applied 
    ! Robotics Kinematics, Dynamics, and Control."

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: dt = 1.0d-3
    integer(int32), parameter :: n = 1000
    real(real64), parameter :: f1 = 5.0d0
    real(real64), parameter :: f2 = 3.0d0
    real(real64), parameter :: f3 = 1.0d1
    real(real64), parameter :: q1 = 0.5d0
    real(real64), parameter :: q2 = 0.25d0
    real(real64), parameter :: q3 = 1.25d0

    ! Linkage Parameters
    real(real64), parameter :: L0 = 1.5d0
    real(real64), parameter :: Alpha0 = -pi / 2.0d0
    real(real64), parameter :: Alpha1 = pi / 2.0d0

    ! Variables
    integer(int32) :: i, jtypes(3)
    real(real64) :: t(n), alpha(n, 3), a(n, 3), theta(n, 3), d(n, 3), X(n, 3), &
        Ti(4, 4), dtheta1(n), dtheta2(n), dd3(n), V(n, 6), J(6, 3), Vfd(n, 3), &
        qdot(3), d3, t1, t2

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: pltX, pltY, pltZ
    type(plot_data_2d) :: pd1, pd2
    class(terminal), pointer :: term
    class(legend), pointer :: lgnd

    ! Define the linkage joint types
    jtypes = [ &
        REVOLUTE_JOINT, &
        REVOLUTE_JOINT, &
        PRISMATIC_JOINT &
    ]

    ! Define the joint motions as functions of time
    t = (/ (dt * i, i = 0, n - 1) /)
    alpha(:,1) = Alpha0
    alpha(:,2) = Alpha1
    alpha(:,3) = 0.0d0
    a = 0.0d0
    theta(:,1) = q1 * sin(2.0d0 * pi * f1 * t)
    theta(:,2) = q2 * sin(2.0d0 * pi * f2 * t)
    theta(:,3) = 0.0d0
    d(:,1:2) = 0.0d0
    d(:,3) = q3 * sin(2.0d0 * pi * f3 * t)

    ! Compute the time derivatives of the joint variables via a finite 
    ! difference
    dtheta1 = finite_difference(t, theta(:,1))
    dtheta2 = finite_difference(t, theta(:,2))
    dd3 = finite_difference(t, d(:,3))

    ! Solve the forward problem and also apply the Jacobian to the velocity
    ! problem
    do i = 1, n
        ! Forward Kinematics Problem
        Ti = dh_forward_kinematics(alpha(i,:), a(i,:), theta(i,:), d(i,:))
        X(i,:) = Ti(1:3,4)

        ! Velocity Problem
        J = dh_jacobian(alpha(i,:), a(i,:), theta(i,:), d(i,:), jtypes)
        qdot = [dtheta1(i), dtheta2(i), dd3(i)]
        V(i,:) = matmul(J, qdot)
    end do

    ! For comparison, compute the time derivatives of the end-effector position
    ! via finite differences
    Vfd(:,1) = finite_difference(t, X(:,1))
    Vfd(:,2) = finite_difference(t, X(:,2))
    Vfd(:,3) = finite_difference(t, X(:,3))

    ! ----------
    ! Plot the results
    call plt%initialize(3, 1)
    call pltX%initialize()
    call pltY%initialize()
    call pltZ%initialize()
    term => plt%get_terminal()
    call term%set_window_height(800)
    call term%set_window_width(1000)
    lgnd => pltX%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call pltX%set_title("x(t)")
    call pltY%set_title("y(t)")
    call pltZ%set_title("z(t)")

    call pd1%define_data(t, V(:,1))
    call pd1%set_name("Analytical")
    call pd1%set_line_width(2.0)
    call pltX%push(pd1)

    call pd2%define_data(t, Vfd(:,1))
    call pd2%set_name("Numerical")
    call pd2%set_line_width(2.0)
    call pd2%set_line_style(LINE_DASHED)
    call pltX%push(pd2)

    ! ----
    call pd1%define_data(t, V(:,2))
    call pltY%push(pd1)

    call pd2%define_data(t, Vfd(:,2))
    call pltY%push(pd2)

    ! ----
    call pd1%define_data(t, V(:,3))
    call pltZ%push(pd1)

    call pd2%define_data(t, Vfd(:,3))
    call pltZ%push(pd2)

    call plt%set(1, 1, pltX)
    call plt%set(2, 1, pltY)
    call plt%set(3, 1, pltZ)
    call plt%draw()
end program