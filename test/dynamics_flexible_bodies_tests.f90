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
module dynamics_flexible_bodies_tests
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_is_finite
    use autodiff
    use dynamics
    implicit none

    type, extends(quadratic_beam_section_energy), private :: coupled_section_energy
        real(real64) :: coupling
    contains
        procedure :: energy_density => coupled_energy_density
    end type

contains
function test_spatial_element_mass() result(ok)
    !! Compares reference mass against closed-form axial, bending and torsional coefficients.
    logical :: ok
    type(flexible_truss_3d) :: truss
    type(flexible_euler_bernoulli_3d) :: eb
    type(flexible_timoshenko_3d) :: tm
    type(flexible_damper_3d) :: damper
    type(flexible_spring_3d) :: spring
    real(real64) :: m6(6,6), m12(12,12)
    real(real64), parameter :: length = 2.0d0, area = 3.0d0, rho = 4.0d0, iy = 0.2d0, iz = 0.4d0
    real(real64), parameter :: mass = rho*area*length, tol = 1.0d-12

    truss%reference_coordinates(1,2) = length
    truss%area = area
    truss%material = linear_elastic(100.0d0,0.3d0,rho)
    m6 = truss%mass_matrix()
    ok = abs(m6(2,2)-mass/3) < tol .and. abs(m6(2,5)-mass/6) < tol
    spring%nodal_mass = [2.0d0, 3.0d0]
    m6 = spring%mass_matrix()
    ok = ok .and. abs(m6(1,1)-2) < tol .and. abs(m6(4,4)-3) < tol .and. abs(m6(1,4)) < tol
    ok = ok .and. maxval(abs(damper%mass_matrix())) < tol

    eb%reference_coordinates = truss%reference_coordinates
    eb%material = truss%material
    eb%area = area
    eb%second_moment_of_area = [iy, iz]
    m12 = eb%mass_matrix()
    ok = ok .and. maxval(abs(m12-transpose(m12))) < tol
    ok = ok .and. abs(m12(2,2) - (mass*13/35 + rho*iz*6/(5*length))) < tol
    ok = ok .and. abs(m12(4,10) - rho*(iy+iz)*length/6) < tol

    tm%reference_coordinates = eb%reference_coordinates
    tm%material = eb%material
    tm%area = area
    tm%second_moment_of_area = eb%second_moment_of_area
    ! Rotate the reference chord to global z; local y stays global y.
    tm%reference_coordinates(:,2) = [0.0d0, 0.0d0, length]
    m12 = tm%mass_matrix()
    ok = ok .and. maxval(abs(m12-transpose(m12))) < tol
    ok = ok .and. abs(m12(6,12) - rho*(iy+iz)*length/6) < tol
    ok = ok .and. abs(m12(4,10) - rho*iz*length/6) < tol .and. abs(m12(2,8)-mass/6) < tol
    if (.not. ok) print *, "TEST FAILED: spatial element mass"
end function

function test_multibody_free_vibration() result(ok)
    !! Checks a fixed floating frame with a single axial truss DOF against
    !! omega=sqrt(3 E/(rho L^2)); the mass uses the element's consistent matrix.
    logical :: ok
    type(multibody_model) :: model
    type(multibody_integrator) :: integrator
    type(multibody_state) :: state
    type(multibody_state), allocatable :: history(:)
    type(flexible_truss_3d) :: truss
    type(rigid_body) :: carrier
    type(dynamic_joint) :: joint
    real(real64) :: q(7), v(7), omega, dt, max_error
    integer :: i

    carrier = rigid_body(m=0.0d0)
    carrier%inertia = 0.0d0
    call model%add_body(carrier, 1)
    truss%reference_coordinates(1,2) = 2.0d0
    truss%material = linear_elastic(12.0d0, density=3.0d0)
    call model%add_element(1, truss, [0,0,0,1,0,0])
    joint = test_joint(FIXED_JOINT, 0, 1)
    call model%add_joint(joint)
    q = 0.0d0
    v = 0.0d0
    q(7) = 1.0d-4
    v(2) = 0.1d0 ! A transverse correction does not change the axial momentum.
    call integrator%initialize(model, state, 0.0d0, q, v)
    omega = sqrt(3.0d0)
    dt = 2*acos(-1.0d0)/(omega*120)
    history = integrator%solve(model, state, dt, 120)
    max_error = 0.0d0
    do i = 1, size(history)
        max_error = max(max_error, abs(history(i)%position(7)/q(7)-cos(omega*history(i)%time)))
    end do
    ok = max_error < 2.0d-3 .and. maxval(abs(history(1)%velocity)) < 1.0d-12
    ok = ok .and. maxval(abs(state%position(:6))) < 1.0d-12 .and. abs(state%time-120*dt) < 1.0d-12
    if (.not. ok) print *, "TEST FAILED: multibody free vibration; normalized error = ", max_error
end function

function test_multibody_elastic_pendulum() result(ok)
    !! A rigid bob is fixed to a deforming spring tip on a revolute carrier.
    !! At the stretched equilibrium, the linear angular mode has frequency
    !! sqrt(m*g*r/(I_carrier+I_bob+m*r^2)).
    logical :: ok
    type(multibody_model) :: model
    type(multibody_integrator) :: integrator
    type(multibody_state) :: state
    type(flexible_spring_3d) :: spring
    type(flexible_damper_3d) :: damper
    type(rigid_body) :: carrier, bob
    type(dynamic_joint) :: joint
    real(real64), parameter :: length = 1.0d0, gravity = 9.81d0, stiffness = 80.0d0, angle = 0.002d0
    real(real64) :: q(13), v(13), radius, omega, dt, max_c, max_angle
    integer :: i

    carrier = rigid_body()
    bob = rigid_body(m=2.0d0)
    call model%add_body(carrier, 1)
    call model%add_body(bob)
    spring%reference_coordinates(2,2) = -length
    spring%stiffness = stiffness
    call model%add_element(1, spring, [0,0,0,0,1,0])
    damper%reference_coordinates = spring%reference_coordinates
    call model%add_element(1, damper, [0,0,0,0,1,0])
    joint = test_joint(REVOLUTE_JOINT, 0, 1)
    call model%add_joint(joint)
    joint = test_joint(FIXED_JOINT, 1, 2)
    joint%parent_frame(2,4) = -length
    call model%add_joint(joint, parent_map=[0,1,0,0,0,0])
    model%gravity = [0.0d0,-gravity,0.0d0]

    radius = length + bob%mass*gravity/stiffness
    omega = sqrt(bob%mass*gravity*radius/(carrier%inertia(3,3)+bob%inertia(3,3)+bob%mass*radius**2))
    q = 0.0d0
    v = 0.0d0
    q(6) = angle
    q(7) = length-radius
    q(8:9) = [radius*sin(angle), -radius*cos(angle)]
    q(13) = angle
    call integrator%initialize(model, state, 0.0d0, q, v)
    dt = 2*acos(-1.0d0)/(omega*180)
    max_c = 0.0d0
    max_angle = 0.0d0
    do i = 1, 180
        call integrator%step(model, state, dt)
        max_c = max(max_c, maxval(abs(model%constraints(state%time,state%position))))
        max_angle = max(max_angle, abs(state%position(6)/angle-cos(omega*state%time)))
    end do
    ok = max_c < 1.0d-10 .and. max_angle < 2.0d-3 .and. abs(state%position(7)-q(7)) < 1.0d-5
    if (.not. ok) print *, "TEST FAILED: elastic pendulum; residual, normalized angular error = ", max_c, max_angle
end function

function test_multibody_prescribed_motion() result(ok)
    !! Exercises time-dependent constraints and projection with an exact accelerating motion.
    logical :: ok
    type(multibody_model) :: model
    type(multibody_integrator) :: integrator
    type(multibody_state) :: state
    type(rigid_body) :: carrier
    real(real64) :: q(6)
    integer :: i

    carrier = rigid_body()
    call model%add_body(carrier)
    call model%set_constraints(1, prescribed_translation)
    q = 0.0d0
    call integrator%initialize(model, state, 0.0d0, q, q)
    ok = abs(state%acceleration(1)-1.0d0) < 1.0d-8
    do i = 1, 10
        call integrator%step(model,state,0.05d0)
    end do
    ok = ok .and. abs(state%position(1)-0.125d0) < 1.0d-10 .and. abs(state%velocity(1)-0.5d0) < 1.0d-8
    ok = ok .and. abs(state%multipliers(1)+carrier%mass) < 1.0d-7
    if (.not. ok) print *, "TEST FAILED: multibody prescribed motion"
end function

subroutine prescribed_translation(time,q,value)
    !! Prescribes x(t)=t^2/2 on the first body.
    real(real64), intent(in) :: time, q(:)
    real(real64), intent(out) :: value(:)
    value(1) = q(1) - 0.5d0*time*time
end subroutine

function test_joint(kind, parent, child) result(joint)
    !! Creates a joint with initially coincident identity attachment frames.
    integer, intent(in) :: kind, parent, child
    type(dynamic_joint) :: joint
    integer :: i

    joint%joint_type = kind
    joint%parent_body = parent
    joint%child_body = child
    joint%parent_frame = 0.0d0
    do i = 1, 4
        joint%parent_frame(i,i) = 1.0d0
    end do
    joint%child_frame = joint%parent_frame
end function

! ------------------------------------------------------------------------------
function test_flexible_truss_2d() result(rst)
    !! Checks the 2D truss against closed-form axial force and stiffness.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-9
    real(real64) :: l, e, fexp
    real(real64), dimension(4) :: u, f
    real(real64), dimension(4,4) :: k, kexp
    type(flexible_truss_2d) :: truss

    rst = .true.
    truss%reference_coordinates(:,1) = [0.0d0, 0.0d0]
    truss%reference_coordinates(:,2) = [1.0d0, 0.0d0]
    truss%area = 2.0d0
    truss%material = linear_elastic(youngs_modulus = 1000.0d0)

    ! Stretched state: L = 1.01, E = (L^2 - L0^2) / (2 L0^2)
    u = [0.0d0, 0.0d0, 0.01d0, 0.0d0]
    call truss%internal_force(u, f, k)
    l = 1.01d0
    e = (l * l - 1.0d0) / 2.0d0
    fexp = 1000.0d0 * 2.0d0 * e * l
    if (abs(f(3) - fexp) > tol .or. abs(f(1) + fexp) > tol .or. &
        abs(f(2)) > tol .or. abs(f(4)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_truss_2d - force"
        rst = .false.
        return
    end if
    if (maxval(abs(k - transpose(k))) > tol) then
        print "(A)", "TEST FAILED: test_flexible_truss_2d - symmetry"
        rst = .false.
        return
    end if

    ! Unstrained state: K = EA/L0 with the axial pattern [1 0 -1 0]
    u = 0.0d0
    call truss%internal_force(u, f, k)
    kexp = 0.0d0
    kexp(1,1) = 2000.0d0
    kexp(3,3) = 2000.0d0
    kexp(1,3) = -2000.0d0
    kexp(3,1) = -2000.0d0
    if (maxval(abs(f)) > tol .or. maxval(abs(k - kexp)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_truss_2d - unstrained"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_rigid_translation() result(rst)
    !! A rigid translation of a 3D truss must produce no internal force.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-10
    real(real64), dimension(6) :: u, f
    real(real64), dimension(6,6) :: k
    type(flexible_truss_3d) :: truss

    rst = .true.
    truss%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    truss%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    truss%area = 1.5d0
    truss%material = linear_elastic(youngs_modulus = 200.0d0)

    u = [1.0d0, 2.0d0, 3.0d0, 1.0d0, 2.0d0, 3.0d0]
    call truss%internal_force(u, f, k)
    if (maxval(abs(f)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_rigid_translation"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_truss_3d_fd() result(rst)
    !! Compares the AD tangent of a nonlinear 3D truss with central differences
    !! of its internal force.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-5
    real(real64), dimension(6) :: u, f
    real(real64), dimension(6,6) :: k, kfd
    type(flexible_truss_3d) :: truss

    rst = .true.
    truss%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    truss%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    truss%area = 1.5d0
    truss%material = polynomial_elastic(linear_elastic(youngs_modulus = 200.0d0), &
        quadratic = 5.0d0, quartic = 2.0d0)

    u = [0.1d0, -0.05d0, 0.2d0, 0.3d0, 0.2d0, -0.1d0]
    call truss%internal_force(u, f, k)
    call fd_stiffness(truss, u, kfd)
    if (maxval(abs(k - kfd)) > tol * max(1.0d0, maxval(abs(k)))) then
        print "(A)", "TEST FAILED: test_flexible_truss_3d_fd"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_spring() result(rst)
    !! Checks a nonlinear 2D spring against closed-form force and stiffness,
    !! then checks a 3D spring against finite differences.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-9
    real(real64) :: delta, fexp, kexp
    real(real64), dimension(4) :: u, f
    real(real64), dimension(4,4) :: k
    real(real64), dimension(6) :: u3, f3
    real(real64), dimension(6,6) :: k3, kfd
    type(flexible_spring_2d) :: spring
    type(flexible_spring_3d) :: spring3

    rst = .true.
    spring%reference_coordinates(:,1) = [0.0d0, 0.0d0]
    spring%reference_coordinates(:,2) = [1.0d0, 0.0d0]
    spring%stiffness = 10.0d0
    spring%cubic_stiffness = 5.0d0

    u = [0.0d0, 0.0d0, 0.2d0, 0.0d0]
    call spring%internal_force(u, f, k)
    delta = 0.2d0
    fexp = 10.0d0 * delta + 5.0d0 * delta**3
    kexp = 10.0d0 + 3.0d0 * 5.0d0 * delta * delta
    if (abs(f(3) - fexp) > tol .or. abs(f(1) + fexp) > tol .or. &
        abs(k(3,3) - kexp) > tol .or. abs(k(1,3) + kexp) > tol) then
        print "(A)", "TEST FAILED: test_flexible_spring - 2D"
        rst = .false.
        return
    end if

    spring3%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    spring3%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    spring3%stiffness = 3.0d0
    spring3%cubic_stiffness = 1.5d0
    u3 = [0.05d0, -0.1d0, 0.02d0, 0.2d0, 0.1d0, -0.3d0]
    call spring3%internal_force(u3, f3, k3)
    call fd_stiffness(spring3, u3, kfd)
    if (maxval(abs(k3 - kfd)) > 1.0d-5 * max(1.0d0, maxval(abs(k3)))) then
        print "(A)", "TEST FAILED: test_flexible_spring - 3D"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_damper() result(rst)
    !! Checks the damping matrix of 2D and 3D viscous dampers at the reference
    !! configuration, where C = c b b^T with b = [-n, n] and n the unit axis.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-10
    real(real64), dimension(4) :: u2
    real(real64), dimension(6) :: u3
    real(real64), dimension(4) :: b2
    real(real64), dimension(6) :: b3
    real(real64), dimension(4,4) :: c2, c2exp
    real(real64), dimension(6,6) :: c3, c3exp
    integer(int32) :: row, col
    type(flexible_damper_2d) :: damper2
    type(flexible_damper_3d) :: damper3

    rst = .true.
    damper2%reference_coordinates(:,1) = [0.0d0, 0.0d0]
    damper2%reference_coordinates(:,2) = [3.0d0, 4.0d0]
    damper2%damping = 2.0d0
    u2 = 0.0d0
    call damper2%damping_matrix(u2, c2)
    b2 = [-0.6d0, -0.8d0, 0.6d0, 0.8d0]
    do col = 1, 4
        do row = 1, 4
            c2exp(row,col) = 2.0d0 * b2(row) * b2(col)
        end do
    end do
    if (maxval(abs(c2 - c2exp)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_damper - 2D"
        rst = .false.
    end if

    damper3%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    damper3%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    damper3%damping = 0.5d0
    u3 = 0.0d0
    call damper3%damping_matrix(u3, c3)
    b3 = [-1.0d0/3.0d0, -2.0d0/3.0d0, -2.0d0/3.0d0, &
        1.0d0/3.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0]
    do col = 1, 6
        do row = 1, 6
            c3exp(row,col) = 0.5d0 * b3(row) * b3(col)
        end do
    end do
    if (maxval(abs(c3 - c3exp)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_damper - 3D"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_beam_linear_limits() result(rst)
    !! Compares unloaded tangents against the classical planar bending and
    !! axial stiffnesses and checks spatial torsional stiffness independently.
    logical :: rst
    type(flexible_euler_bernoulli_2d) :: eb
    type(flexible_timoshenko_2d) :: tm
    type(flexible_euler_bernoulli_3d) :: eb3
    type(flexible_timoshenko_3d) :: tm3
    real(real64) :: f(6), k(6,6), f3(12), k3(12,12), expected(2,2)
    integer, parameter :: tip(2) = [5, 6]
    real(real64), parameter :: l = 2.0d0, ei = 3.0d0, ga = 7.0d0, gj = 11.0d0
    real(real64), parameter :: e = 5.0d0, g = 2.5d0
    real(real64), parameter :: tol = 1.0d-10, twist = 0.6d0
    real(real64) :: u3(12)

    eb%reference_coordinates(:,2) = [l, 0.0d0]
    eb%second_moment_of_area = ei/e
    eb%area = 2.0d0
    eb%material = linear_elastic(youngs_modulus = e, poissons_ratio = 0.0d0)
    tm%reference_coordinates = eb%reference_coordinates
    tm%second_moment_of_area = eb%second_moment_of_area
    tm%effective_shear_area = ga/g
    tm%area = eb%area
    tm%material = eb%material

    call eb%internal_force([0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], f, k)
    expected = ei * reshape([12/l**3, -6/l**2, -6/l**2, 4/l], [2,2])
    rst = all(ieee_is_finite(k)) .and. maxval(abs(k(tip,tip) - expected)) < tol .and. &
        abs(k(4,4) - 10/l) < tol .and. maxval(abs(f)) < tol
    if (.not. rst) then
        print *, "TEST FAILED: Euler-Bernoulli planar linear stiffness"
        return
    end if

    call tm%internal_force([0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], f, k)
    expected = reshape([ga/l, -ga/2, -ga/2, ei/l + ga*l/4], [2,2])
    rst = all(ieee_is_finite(k)) .and. maxval(abs(k(tip,tip) - expected)) < tol .and. &
        abs(k(4,4) - 10/l) < tol .and. maxval(abs(f)) < tol
    if (.not. rst) then
        print *, "TEST FAILED: Timoshenko planar linear stiffness"
        return
    end if

    eb3%reference_coordinates(:,2) = [l, 0.0d0, 0.0d0]
    eb3%torsion_constant = gj/g
    eb3%material = eb%material
    tm3%reference_coordinates = eb3%reference_coordinates
    tm3%torsion_constant = eb3%torsion_constant
    tm3%material = eb%material
    u3 = 0.0d0
    call eb3%internal_force(u3, f3, k3)
    rst = all(ieee_is_finite(k3)) .and. abs(k3(10,10) - gj/l) < tol
    u3(10) = twist
    call tm3%internal_force(u3, f3, k3)
    rst = rst .and. all(ieee_is_finite(k3)) .and. abs(f3(10) - gj*twist/l) < tol .and. &
        abs(k3(10,10) - gj/l) < tol
    if (.not. rst) print *, "TEST FAILED: spatial beam torsion"
    if (rst) rst = test_flexible_beam_material_swap()
end function

function test_flexible_beam_material_swap() result(rst)
    logical :: rst
    type(flexible_timoshenko_3d) :: beam
    real(real64), parameter :: l = 2.0d0, tol = 1.0d-10
    real(real64), parameter :: e(2) = [12.0d0, 18.0d0], g(2) = [5.0d0, 8.0d0]
    real(real64), parameter :: nu(2) = [0.2d0, 0.125d0]
    real(real64) :: u(12), f(12), k(12,12), expected(6)
    integer :: i, j

    beam%reference_coordinates(:,2) = [l, 0.0d0, 0.0d0]
    beam%area = 2.0d0
    beam%second_moment_of_area = [3.0d0, 4.0d0]
    beam%torsion_constant = 5.0d0
    beam%effective_shear_area = [0.5d0, 1.5d0]
    u = 0.0d0
    do i = 1, 2
        beam%material = linear_elastic(youngs_modulus = e(i), poissons_ratio = nu(i))
        call beam%internal_force(u, f, k)
        ! Tip rotation about y uses shear along z, and vice versa.
        expected = [2*e(i)/l, 0.5d0*g(i)/l, 1.5d0*g(i)/l, 5*g(i)/l, &
            3*e(i)/l + 1.5d0*g(i)*l/4, 4*e(i)/l + 0.5d0*g(i)*l/4]
        rst = all(ieee_is_finite(k)) .and. maxval(abs([(k(j,j), j=7,12)] - expected)) < tol
        if (.not. rst) then
            print *, "TEST FAILED: spatial beam material swap", i
            return
        end if
    end do
end function

function coupled_energy_density(this, geometry, material, strain) result(w)
    class(coupled_section_energy), intent(in) :: this
    type(beam_section_geometry), intent(in) :: geometry
    class(flexible_material), intent(in) :: material
    type(beam_section_strain), intent(in) :: strain
    type(hyperdual) :: w

    w = this%quadratic_beam_section_energy%energy_density(geometry, material, strain) + &
        this%coupling * geometry%area * strain%extension * strain%curvature(1)
end function

function test_flexible_beam_section_energy() result(rst)
    logical :: rst
    type(flexible_euler_bernoulli_3d) :: eb
    type(flexible_timoshenko_3d) :: tm
    real(real64), parameter :: l = 2.0d0, area = 1.5d0, coupling = 3.0d0, delta = 0.1d0
    real(real64), parameter :: cross_tangent = area * coupling / l, tol = 1.0d-10
    real(real64) :: u(12), f(12), k(12,12)

    eb%reference_coordinates(1,2) = l
    eb%area = area
    eb%material = linear_elastic(youngs_modulus = 20.0d0)
    allocate(eb%section_law, source = coupled_section_energy(coupling = coupling))
    tm%reference_coordinates = eb%reference_coordinates
    tm%area = eb%area
    tm%material = eb%material
    allocate(tm%section_law, source = eb%section_law)
    u = 0.0d0
    u(7) = delta

    call eb%internal_force(u, f, k)
    rst = abs(f(10) - cross_tangent*delta) < tol .and. abs(k(10,7) - cross_tangent) < tol
    if (.not. rst) print *, "TEST FAILED: Euler-Bernoulli section energy coupling"

    call tm%internal_force(u, f, k)
    if (.not. (abs(f(10) - cross_tangent*delta) < tol .and. abs(k(10,7) - cross_tangent) < tol)) then
        print *, "TEST FAILED: Timoshenko section energy coupling"
        rst = .false.
    end if
end function

function test_flexible_beam_objectivity() result(rst)
    !! Verifies that a finite rigid motion gives no elastic force and that
    !! superposing that motion on a strained state preserves its energy.
    logical :: rst
    type(flexible_euler_bernoulli_3d) :: eb
    type(flexible_timoshenko_3d) :: tm
    real(real64), parameter :: axis(3) = [2.0d0, -1.0d0, 2.0d0]/3.0d0
    real(real64), parameter :: shift(3) = [0.2d0, -0.5d0, 0.4d0], theta = 1.1d0
    real(real64) :: u(12), rigid(12), moved(12), p(3), r(3), unit_axis(3), angle
    type(quaternion) :: superposed
    integer :: i, j

    eb%reference_coordinates(:,1) = [0.2d0, 0.1d0, -0.4d0]
    eb%reference_coordinates(:,2) = [1.5d0, 1.0d0, 0.6d0]
    eb%second_moment_of_area = [2.0d0, 5.0d0]/20.0d0
    eb%torsion_constant = 1.0d0/8.0d0
    eb%material = polynomial_elastic(linear_elastic(20.0d0, 0.25d0), quadratic = 9.0d0, quartic = 3.0d0)
    tm%reference_coordinates = eb%reference_coordinates
    tm%second_moment_of_area = eb%second_moment_of_area
    tm%torsion_constant = eb%torsion_constant
    tm%effective_shear_area = [3.0d0, 7.0d0]/8.0d0
    tm%material = eb%material
    u = [0.01d0, 0.04d0, -0.03d0, -0.15d0, 0.1d0, 0.05d0, &
        0.2d0, -0.08d0, 0.1d0, 0.2d0, 0.15d0, -0.25d0]
    do i = 1, 2
        j = 6*(i-1)
        p = eb%reference_coordinates(:,i)
        rigid(j+1:j+3) = rotate_test_point(p, axis, theta) + shift - p
        rigid(j+4:j+6) = theta*axis
        moved(j+1:j+3) = rotate_test_point(p + u(j+1:j+3), axis, theta) + shift - p
        r = u(j+4:j+6)
        superposed = quaternion(theta, axis) * quaternion(norm2(r), r/norm2(r))
        call superposed%to_angle_axis(angle, unit_axis)
        moved(j+4:j+6) = angle*unit_axis
    end do
    rst = check_beam_objectivity(eb, u, rigid, moved) .and. check_beam_objectivity(tm, u, rigid, moved)
    if (.not. rst) print *, "TEST FAILED: flexible beam finite-motion objectivity"
end function

function rotate_test_point(p, axis, angle) result(q)
    !! Uses Rodrigues' formula to construct the independent rigid-motion test.
    real(real64), intent(in) :: p(3), axis(3), angle
    real(real64) :: q(3)

    q = cos(angle)*p + sin(angle)*cross_product(axis, p) + (1-cos(angle))*dot_product(axis,p)*axis
end function

function check_beam_objectivity(el, u, rigid, moved) result(rst)
    !! Checks energy invariance and compares the tangent at a strained state
    !! with central differences of the internal force.
    class(finite_element), intent(in) :: el
    real(real64), intent(in) :: u(12), rigid(12), moved(12)
    logical :: rst
    real(real64) :: f(12), k(12,12), kfd(12,12)
    type(hyperdual) :: seed(12), w, moved_w, rigid_w
    integer :: i

    do i = 1, 12
        seed(i) = hyperdual(u(i), 0.0d0, 0.0d0, 0.0d0)
    end do
    w = el%strain_energy(seed)
    seed%x = moved
    moved_w = el%strain_energy(seed)
    seed%x = rigid
    rigid_w = el%strain_energy(seed)
    call el%internal_force(rigid, f, k)
    rst = el%get_dof_count() == 12 .and. abs(rigid_w%x) < 1.0d-12 .and. &
        maxval(abs(f)) < 1.0d-10 .and. w%x > 1.0d-5 .and. abs(w%x - moved_w%x) < 1.0d-11
    call el%internal_force(moved, f, k)
    call fd_stiffness(el, moved, kfd)
    rst = rst .and. all(ieee_is_finite(k)) .and. all(ieee_is_finite(kfd)) .and. &
        maxval(abs(k - kfd)) < 1.0d-6 * max(1.0d0, maxval(abs(k)))
end function

! ------------------------------------------------------------------------------
subroutine fd_stiffness(el, u, kfd)
    !! Estimates the tangent stiffness by central differences of the internal
    !! force of a finite element.
    class(finite_element), intent(in) :: el
    real(real64), intent(in), dimension(:) :: u
    real(real64), intent(out), dimension(:,:) :: kfd

    real(real64), parameter :: h = 1.0d-6
    integer(int32) :: j, n
    real(real64), allocatable, dimension(:) :: up, fp, fm
    real(real64), allocatable, dimension(:,:) :: kdummy

    n = size(u)
    allocate(up(n), fp(n), fm(n), kdummy(n,n))
    do j = 1, n
        up = u
        up(j) = u(j) + h
        call el%internal_force(up, fp, kdummy)
        up = u
        up(j) = u(j) - h
        call el%internal_force(up, fm, kdummy)
        kfd(:,j) = (fp - fm) / (2.0d0 * h)
    end do
end subroutine

! ------------------------------------------------------------------------------
end module
