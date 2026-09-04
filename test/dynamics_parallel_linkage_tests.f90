module dynamics_parallel_linkage_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_graph
    use dynamics_joints
    use dynamics_linkage
    use dynamics_parallel_linkage
    use dynamics_kinematics
    implicit none

    ! Four-bar linkage geometry (a Grashof crank-rocker)
    real(real64), parameter :: fb_crank = 1.0d0
    real(real64), parameter :: fb_coupler = 3.5d0
    real(real64), parameter :: fb_rocker = 3.0d0
    real(real64), parameter :: fb_ground = 4.0d0

contains
! ------------------------------------------------------------------------------
    pure function translate_x(a) result(rst)
        ! A 4-by-4 translation along the x-axis.
        real(real64), intent(in) :: a
        real(real64) :: rst(4, 4)
        integer(int32) :: i
        rst = 0.0d0
        do i = 1, 4
            rst(i,i) = 1.0d0
        end do
        rst(1,4) = a
    end function

! ------------------------------------------------------------------------------
    function planar_link(length) result(rst)
        ! Constructs a planar link of the requested length carrying a joint
        ! frame at each end.  The body frame is coincident with the first joint.
        real(real64), intent(in) :: length
        type(multi_joint_link) :: rst
        real(real64) :: frames(4, 4, 2)
        frames(:,:,1) = translate_x(0.0d0)
        frames(:,:,2) = translate_x(length)
        rst = multi_joint_link(frames)
    end function

! ------------------------------------------------------------------------------
    function build_four_bar() result(rst)
        !! Builds a planar four-bar linkage.  The ground pivots lie at the
        !! origin and at (fb_ground, 0); the crank is driven about the origin.
        !! The end-effector is the distal end of the coupler.
        type(planar_linkage) :: rst

        ! Local Variables
        type(link_container) :: lnks(4)
        type(joint) :: jnts(4)

        ! Define the links
        allocate(lnks(1)%item, source = planar_link(fb_ground))
        allocate(lnks(2)%item, source = planar_link(fb_crank))
        allocate(lnks(3)%item, source = planar_link(fb_coupler))
        allocate(lnks(4)%item, source = planar_link(fb_rocker))

        ! Define the joints
        jnts(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, .true.)
        jnts(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
        jnts(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
        jnts(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)

        ! Build the mechanism
        rst = planar_linkage(lnks, jnts, base = 1, effector = 3, &
            tool = translate_x(fb_coupler))
    end function

! ------------------------------------------------------------------------------
    subroutine four_bar_solution(theta, q, c)
        !! Computes the closed-form solution of the four-bar linkage.
        real(real64), intent(in) :: theta
            !! The crank angle.
        real(real64), intent(out) :: q(4)
            !! The corresponding joint variables.
        real(real64), intent(out) :: c(2)
            !! The location of the coupler-rocker pivot.

        ! Local Variables
        real(real64) :: b(2), d(2), u(2), w(2), L, x, h, psi, phi

        ! Locate the crank tip and the fixed pivot
        b = [fb_crank * cos(theta), fb_crank * sin(theta)]
        d = [fb_ground, 0.0d0]

        ! Intersect the coupler and rocker circles
        L = norm2(d - b)
        u = (d - b) / L
        w = [-u(2), u(1)]
        x = 0.5d0 * (L**2 + fb_coupler**2 - fb_rocker**2) / L
        h = sqrt(fb_coupler**2 - x**2)
        c = b + x * u + h * w

        ! Assemble the joint variables
        psi = atan2(c(2) - b(2), c(1) - b(1))
        phi = atan2(c(2) - d(2), c(1) - d(1))
        q = [theta, psi - theta, phi - psi, -phi]
    end subroutine

! ------------------------------------------------------------------------------
    function test_graph_topology() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        type(graph) :: g
        type(spanning_tree) :: tree
        type(graph_loop), allocatable, dimension(:) :: loops
        type(graph_path) :: path
        integer(int32), allocatable, dimension(:) :: cuts

        ! Initialization
        rst = .true.

        ! Build the topology of a four-bar linkage
        call g%initialize(4)
        call g%add_edge(1, 2)
        call g%add_edge(2, 3)
        call g%add_edge(3, 4)
        call g%add_edge(4, 1)

        ! Test the basic counts
        if (.not.assert(g%get_vertex_count(), 4)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -1"
        end if
        if (.not.assert(g%get_edge_count(), 4)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -2"
        end if
        if (.not.assert(g%get_independent_loop_count(), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -3"
        end if

        ! Test the spanning tree
        tree = g%build_spanning_tree(1)
        if (.not.tree%is_connected()) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -4"
        end if
        cuts = tree%get_cut_edges()
        if (.not.assert(size(cuts), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -5"
        end if
        loops = g%find_independent_loops(tree)
        if (.not.assert(size(loops), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -6"
        end if

        ! Every vertex is at most two edges from the root
        path = tree%get_path(3)
        if (.not.assert(size(path%edges), 2)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -7"
        end if
        path = tree%get_path(1)
        if (.not.assert(size(path%edges), 0)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_graph_topology -8"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_four_bar_topology() result(rst)
        ! Arguments
        logical :: rst

        ! Local Variables
        type(planar_linkage) :: linkage

        ! Initialization
        rst = .true.
        linkage = build_four_bar()

        ! Test
        if (.not.assert(linkage%get_link_count(), 4)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -1"
        end if
        if (.not.assert(linkage%get_joint_count(), 4)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -2"
        end if
        if (.not.assert(linkage%get_variable_count(), 4)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -3"
        end if
        if (.not.assert(linkage%get_loop_count(), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -4"
        end if
        if (.not.assert(linkage%get_constraint_count(), 3)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -5"
        end if
        if (.not.assert(linkage%get_degrees_of_freedom(), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -6"
        end if
        if (.not.assert(linkage%get_actuated_variable_count(), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_topology -7"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_four_bar_constraints() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-10

        ! Local Variables
        type(planar_linkage) :: linkage
        real(real64) :: theta, q(4), c(2)
        real(real64), allocatable, dimension(:) :: f, zeros

        ! Initialization
        rst = .true.
        linkage = build_four_bar()
        theta = 0.6d0
        call four_bar_solution(theta, q, c)

        ! The closed-form solution must satisfy the loop-closure constraints
        f = linkage%constraints(q)
        allocate(zeros(size(f)), source = 0.0d0)
        if (.not.assert(f, zeros, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_constraints -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_four_bar_forward_kinematics() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        type(planar_linkage) :: linkage
        real(real64) :: theta, q(4), c(2), p(3), T(4, 4)
        real(real64), allocatable, dimension(:) :: f

        ! Initialization
        rst = .true.
        linkage = build_four_bar()

        ! Establish the assembly mode, then allow the solver to track the
        ! mechanism as the crank is swept
        call four_bar_solution(0.7d0, q, c)
        call linkage%set_configuration(q)

        ! Sweep the crank and compare against the closed-form solution
        do i = 1, 8
            theta = 0.7d0 * real(i, real64)
            call four_bar_solution(theta, q, c)
            T = linkage%forward_kinematics([theta])
            if (.not.assert(T(1:2,4), c, tol)) then
                rst = .false.
                print "(A)", "TEST FAILED: test_four_bar_forward_kinematics -1"
                exit
            end if

            ! The solution must also satisfy the constraints
            f = linkage%constraints(linkage%get_configuration())
            if (maxval(abs(f)) > tol) then
                rst = .false.
                print "(A)", "TEST FAILED: test_four_bar_forward_kinematics -2"
                exit
            end if

            ! The planar pose accessor must agree with the transformation matrix
            p = linkage%end_effector_pose([theta])
            if (.not.assert(p(1:2), c, tol)) then
                rst = .false.
                print "(A)", "TEST FAILED: test_four_bar_forward_kinematics -3"
                exit
            end if
        end do
    end function

! ------------------------------------------------------------------------------
    function test_four_bar_jacobian() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-5
        real(real64), parameter :: h = 1.0d-5

        ! Local Variables
        type(planar_linkage) :: linkage
        real(real64) :: theta, q(4), c(2), p1(3), p2(3), fd(3), jw(3), s, cs
        real(real64), allocatable, dimension(:,:) :: J

        ! Initialization
        rst = .true.
        linkage = build_four_bar()
        theta = 0.85d0
        call four_bar_solution(theta, q, c)
        call linkage%set_configuration(q)

        ! Compute the Jacobian
        J = linkage%jacobian([theta])
        if (.not.assert(size(J, 1), 3) .or. .not.assert(size(J, 2), 1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_jacobian -1"
            return
        end if

        ! The Jacobian is expressed in the end-effector frame; rotate the
        ! translational terms into the frame of the base link
        p1 = linkage%end_effector_pose([theta])
        s = sin(p1(3))
        cs = cos(p1(3))
        jw(1) = cs * J(1,1) - s * J(2,1)
        jw(2) = s * J(1,1) + cs * J(2,1)
        jw(3) = J(3,1)

        ! Compare against a finite difference estimate of the end-effector pose
        call four_bar_solution(theta + h, q, c)
        call linkage%set_configuration(q)
        p1 = linkage%end_effector_pose([theta + h])
        call four_bar_solution(theta - h, q, c)
        call linkage%set_configuration(q)
        p2 = linkage%end_effector_pose([theta - h])
        fd = (p1 - p2) / (2.0d0 * h)

        if (.not.assert(jw, fd, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_jacobian -2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_four_bar_inverse_kinematics() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-6

        ! Local Variables
        type(planar_linkage) :: linkage
        real(real64) :: theta, q(4), c(2), T(4, 4)
        real(real64), allocatable, dimension(:) :: qa

        ! Initialization
        rst = .true.
        linkage = build_four_bar()
        theta = 1.1d0
        call four_bar_solution(theta, q, c)
        call linkage%set_configuration(q)
        T = linkage%forward_kinematics([theta])

        ! Perturb the configuration and solve for the crank angle
        call four_bar_solution(theta - 0.1d0, q, c)
        call linkage%set_configuration(q)
        qa = linkage%inverse_kinematics(T)

        if (.not.assert(qa(1), theta, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_four_bar_inverse_kinematics -1"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_open_chain_equivalence() result(rst)
        !! Verifies that an open-chain mechanism described by a parallel_linkage
        !! reproduces the results of the equivalent serial_linkage.
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        integer(int32) :: i
        type(binary_link) :: links(3)
        type(link_container) :: lnks(4)
        type(joint) :: jnts(3)
        type(serial_linkage) :: serial
        type(parallel_linkage) :: mechanism
        real(real64) :: q(3), Ts(4, 4), Tp(4, 4)

        ! Initialization
        rst = .true.
        links(1) = binary_link(twist = -0.5d0 * pi, offset = 1.5d0)
        links(2) = binary_link(twist = 0.5d0 * pi, length = 0.75d0)
        links(3) = binary_link(offset = 1.25d0, jtype = PRISMATIC_JOINT)
        serial = serial_linkage(links)

        ! The equivalent mechanism requires a ground link
        allocate(lnks(1)%item, source = binary_link())
        do i = 1, 3
            allocate(lnks(i+1)%item, source = links(i))
            jnts(i) = joint(links(i)%joint_type, i, i + 1, 2, 1, .true.)
        end do
        mechanism = parallel_linkage(lnks, jnts, base = 1, effector = 4)

        ! An open chain has no loops
        if (.not.assert(mechanism%get_loop_count(), 0)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_open_chain_equivalence -1"
        end if
        if (.not.assert(mechanism%get_degrees_of_freedom(), 3)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_open_chain_equivalence -2"
        end if

        ! Compare the forward kinematics
        call random_number(q)
        Ts = serial%forward_kinematics(q)
        Tp = mechanism%forward_kinematics(q)
        if (.not.assert(Tp, Ts, tol)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_open_chain_equivalence -3"
        end if
    end function

! ------------------------------------------------------------------------------
end module
