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
module dynamics_flexible_beams
    !! Straight, two-node reference beams with objective finite-rotation
    !! kinematics. All section stiffnesses refer to the reference configuration.
    use iso_fortran_env, only : real64, int32
    use autodiff
    use dynamics_helper, only : cross_product
    use dynamics_flexible_types, only : finite_element
    use dynamics_flexible_materials, only : flexible_material
    implicit none
    private
    public :: flexible_euler_bernoulli_2d, flexible_euler_bernoulli_3d
    public :: flexible_timoshenko_2d, flexible_timoshenko_3d

    real(real64), parameter :: frame_tol = 100.0d0 * epsilon(1.0d0)

    type, abstract, extends(finite_element) :: beam_section
        real(real64) :: area = 1.0d0
            !! Reference cross-sectional area.
        class(flexible_material), allocatable :: material
            !! Required axial law, evaluated at the generalized axial strain
            !! `stretch - 1` rather than the truss Green-Lagrange strain.
    end type

    type, abstract, extends(beam_section) :: planar_beam
        real(real64) :: reference_coordinates(2,2) = 0.0d0
            !! One column per node. DOFs are `[ux, uy, theta_z]` per node,
            !! with positive angles measured counterclockwise in radians.
        real(real64) :: bending_stiffness = 1.0d0
            !! Section stiffness EI about the out-of-plane axis.
    contains
        procedure :: get_dof_count => planar_beam_dof_count
    end type

    type, abstract, extends(beam_section) :: spatial_beam
        real(real64) :: reference_coordinates(3,2) = 0.0d0
            !! One column per node. DOFs are `[ux, uy, uz, rx, ry, rz]` per node.
            !! The last three are a global rotation vector in radians, applied
            !! through the exponential map to the reference section frame.
        real(real64) :: reference_y(3) = [0.0d0, 1.0d0, 0.0d0]
            !! Direction projected onto the reference cross section to set its
            !! local y axis; it must not be parallel to the reference chord.
        real(real64) :: bending_stiffness(2) = 1.0d0
            !! Section stiffnesses [EI_y, EI_z].
        real(real64) :: torsional_stiffness = 1.0d0
            !! Section stiffness GJ about the beam axis.
    contains
        procedure :: get_dof_count => spatial_beam_dof_count
    end type

    type, extends(planar_beam) :: flexible_euler_bernoulli_2d
        !! Shear-free planar beam with a cubic centerline.
    contains
        procedure :: strain_energy => euler_bernoulli_2d_energy
    end type

    type, extends(spatial_beam) :: flexible_euler_bernoulli_3d
        !! Shear-free spatial beam with a cubic centerline and section twist.
    contains
        procedure :: strain_energy => euler_bernoulli_3d_energy
    end type

    type, extends(planar_beam) :: flexible_timoshenko_2d
        !! Extensible planar beam with independent section rotations and shear.
        real(real64) :: shear_stiffness = 1.0d0
            !! Effective section shear stiffness, including the shear correction factor.
    contains
        procedure :: strain_energy => timoshenko_2d_energy
    end type

    type, extends(spatial_beam) :: flexible_timoshenko_3d
        !! Extensible spatial beam with independent section rotations, shear and twist.
        real(real64) :: shear_stiffness(2) = 1.0d0
            !! Effective shear stiffnesses along local [y, z].
    contains
        procedure :: strain_energy => timoshenko_3d_energy
    end type

contains

function planar_beam_dof_count(this) result(n)
    !! Returns three DOFs per node for a two-node planar beam.
    class(planar_beam), intent(in) :: this
        !! The planar element.
    integer(int32) :: n
        !! Total number of active DOFs.
    n = 6
end function

function spatial_beam_dof_count(this) result(n)
    !! Returns six DOFs per node for a two-node spatial beam.
    class(spatial_beam), intent(in) :: this
        !! The spatial element.
    integer(int32) :: n
        !! Total number of active DOFs.
    n = 12
end function

function euler_bernoulli_2d_energy(this, u) result(w)
    !! Evaluates the planar shear-free energy by embedding the motion in 3D.
    class(flexible_euler_bernoulli_2d), intent(in) :: this
    type(hyperdual), intent(in) :: u(:)
        !! Six nodal translations and angles in node order.
    type(hyperdual) :: w
        !! Total strain energy.
    type(hyperdual) :: u3(12)
    real(real64) :: x(3,2), y(3)

    call embed_planar(this%reference_coordinates, u, x, y, u3)
    w = euler_bernoulli_energy(this, x, y, [0.0d0, 0.0d0, this%bending_stiffness], u3)
end function

function euler_bernoulli_3d_energy(this, u) result(w)
    !! Evaluates the spatial axial, bending and torsional energy.
    class(flexible_euler_bernoulli_3d), intent(in) :: this
    type(hyperdual), intent(in) :: u(:)
        !! Twelve nodal translations and rotation-vector components in node order.
    type(hyperdual) :: w
        !! Total strain energy.

    w = euler_bernoulli_energy(this, this%reference_coordinates, this%reference_y, &
        [this%torsional_stiffness, this%bending_stiffness], u)
end function

function timoshenko_2d_energy(this, u) result(w)
    !! Evaluates the planar axial, bending and shear energy.
    class(flexible_timoshenko_2d), intent(in) :: this
    type(hyperdual), intent(in) :: u(:)
        !! Six nodal translations and angles in node order.
    type(hyperdual) :: w
        !! Total strain energy.
    type(hyperdual) :: u3(12)
    real(real64) :: x(3,2), y(3)

    call embed_planar(this%reference_coordinates, u, x, y, u3)
    w = timoshenko_energy(this, x, y, [0.0d0, 0.0d0, this%bending_stiffness], &
        [this%shear_stiffness, 0.0d0], u3)
end function

function timoshenko_3d_energy(this, u) result(w)
    !! Evaluates the spatial axial, bending, shear and torsional energy.
    class(flexible_timoshenko_3d), intent(in) :: this
    type(hyperdual), intent(in) :: u(:)
        !! Twelve nodal translations and rotation-vector components in node order.
    type(hyperdual) :: w
        !! Total strain energy.

    w = timoshenko_energy(this, this%reference_coordinates, this%reference_y, &
        [this%torsional_stiffness, this%bending_stiffness], this%shear_stiffness, u)
end function

subroutine embed_planar(x, u, x3, y, u3)
    !! Embeds the two planar reference nodes and six active DOFs in the xy
    !! plane of the spatial formulation. The reference local y axis is left
    !! of the chord so local z always coincides with global z.
    real(real64), intent(in) :: x(2,2)
        !! Planar reference nodes, one per column.
    type(hyperdual), intent(in) :: u(:)
        !! Six planar DOFs.
    real(real64), intent(out) :: x3(3,2)
        !! Reference nodes embedded in the global xy plane.
    real(real64), intent(out) :: y(3)
        !! Reference section orientation hint.
    type(hyperdual), intent(out) :: u3(12)
        !! Spatial DOFs; the six inactive components are zero.
    integer :: map(6), i

    if (size(u) /= 6) error stop "Planar flexible beam requires six DOFs"
    x3(1:2,:) = x
    x3(3,:) = 0.0d0
    y = [x(2,1) - x(2,2), x(1,2) - x(1,1), 0.0d0]
    map = [1, 2, 6, 7, 8, 12]
    do i = 1, 12
        u3(i) = 0.0d0
    end do
    do i = 1, 6
        u3(map(i)) = u(i)
    end do
end subroutine

subroutine beam_kinematics(section, x, y, u, l0, frame, chord, q)
    !! Validates the section and reference geometry, then returns the reference
    !! length/frame, current chord, and the two global rotation quaternions.
    class(beam_section), intent(in) :: section
    real(real64), intent(in) :: x(3,2), y(3)
        !! Reference coordinates and local-y hint.
    type(hyperdual), intent(in) :: u(:)
        !! Spatial nodal DOFs, including hyper-dual rotation components.
    real(real64), intent(out) :: l0, frame(3,3)
        !! Reference length and orthonormal directors (columns).
    type(hyperdual), intent(out) :: chord(3), q(4,2)
        !! Current end-to-end vector and scalar-first nodal quaternions.
    real(real64) :: dy(3), yn
    integer :: i

    if (size(u) /= 12) error stop "Spatial flexible beam requires twelve DOFs"
    if (.not. allocated(section%material)) error stop "Flexible beam requires an axial material"
    if (.not. (section%area > 0.0d0)) error stop "Flexible beam requires positive area"
    frame(:,1) = x(:,2) - x(:,1)
    l0 = norm2(frame(:,1))
    if (.not. (l0 > 0.0d0)) error stop "Flexible beam has zero reference length"
    frame(:,1) = frame(:,1) / l0
    yn = norm2(y)
    if (.not. (yn > 0.0d0)) error stop "Flexible beam requires a nonzero reference y direction"
    dy = y / yn - dot_product(y / yn, frame(:,1)) * frame(:,1)
    yn = norm2(dy)
    if (.not. (yn > frame_tol)) error stop "Flexible beam reference y is parallel to its chord"
    frame(:,2) = dy / yn
    frame(:,3) = cross_product(frame(:,1), frame(:,2))
    do i = 1, 3
        chord(i) = (x(i,2) - x(i,1)) + (u(6+i) - u(i))
    end do
    q(:,1) = rotation_quaternion(u(4:6))
    q(:,2) = rotation_quaternion(u(10:12))
end subroutine

function euler_bernoulli_energy(section, x, y, stiffness, u) result(w)
    !! Integrates the shear-free energy with four Gauss points on the reference
    !! interval. The nodal section axes set the Hermite tangent directions;
    !! the current chord length scales both end derivatives. The local y axis
    !! is interpolated then projected normal to the current centerline.
    class(beam_section), intent(in) :: section
    real(real64), intent(in) :: x(3,2), y(3), stiffness(3)
        !! Reference geometry, local-y hint, and [GJ, EI_y, EI_z].
    type(hyperdual), intent(in) :: u(:)
        !! Spatial nodal DOFs.
    type(hyperdual) :: w
        !! Total axial, bending and torsional strain energy.
    real(real64), parameter :: point(4) = [0.06943184420297371d0, 0.33000947820757187d0, &
        0.66999052179242813d0, 0.93056815579702629d0]
    real(real64), parameter :: weight(4) = [0.17392742256872693d0, 0.32607257743127307d0, &
        0.32607257743127307d0, 0.17392742256872693d0]
    real(real64) :: l0, frame(3,3), s
    type(hyperdual) :: chord(3), q(4,2), a(3,2), b(3,2), scale, speed, projection, width
    type(hyperdual) :: dr(3), ddr(3), t(3), dt(3), guide(3), dguide(3), d2(3), d3(3), k(3), density
    integer :: i, j

    if (.not. all(stiffness >= 0.0d0)) error stop "Flexible beam stiffness cannot be negative"
    call beam_kinematics(section, x, y, u, l0, frame, chord, q)
    scale = hd_dot(chord, chord)
    if (.not. (scale%x > 0.0d0)) error stop "Euler-Bernoulli beam has zero current chord"
    scale = sqrt(scale)
    do j = 1, 2
        a(:,j) = rotate_axis(q(:,j), frame(:,1))
        b(:,j) = rotate_axis(q(:,j), frame(:,2))
    end do
    do i = 1, 3
        dguide(i) = (b(i,2) - b(i,1)) / l0
    end do

    w = 0.0d0
    do j = 1, 4
        s = point(j)
        do i = 1, 3
            dr(i) = (6*s*(1-s)*chord(i) + scale*((1-4*s+3*s*s)*a(i,1) + (3*s*s-2*s)*a(i,2))) / l0
            ddr(i) = (6*(1-2*s)*chord(i) + scale*((6*s-4)*a(i,1) + (6*s-2)*a(i,2))) / l0**2
            guide(i) = (1-s)*b(i,1) + s*b(i,2)
        end do
        speed = hd_dot(dr, dr)
        if (.not. (speed%x > frame_tol**2)) error stop "Euler-Bernoulli beam has a singular tangent"
        speed = sqrt(speed)
        do i = 1, 3
            t(i) = dr(i) / speed
        end do
        projection = hd_dot(t, ddr)
        do i = 1, 3
            dt(i) = (ddr(i) - projection*t(i)) / speed
        end do
        projection = hd_dot(t, guide)
        do i = 1, 3
            d2(i) = guide(i) - projection*t(i)
        end do
        width = hd_dot(d2, d2)
        if (.not. (width%x > frame_tol**2)) error stop "Euler-Bernoulli beam has a singular section frame"
        width = sqrt(width)
        do i = 1, 3
            d2(i) = d2(i) / width
        end do
        d3 = hd_cross(t, d2)
        ! d3 dot d2' eliminates terms parallel to t and d2.
        k(1) = (hd_dot(d3, dguide) - projection*hd_dot(d3, dt)) / width
        k(2) = -hd_dot(d3, dt)
        k(3) = hd_dot(d2, dt)
        density = section%area * section%material%strain_energy_density(speed - 1.0d0)
        do i = 1, 3
            density = density + 0.5d0 * stiffness(i) * k(i)*k(i)
        end do
        w = w + weight(j) * l0 * density
    end do
end function

function timoshenko_energy(section, x, y, stiffness, shear, u) result(w)
    !! Uses a linear centerline and shortest-arc section interpolation. The
    !! relative rotation gives constant material curvature; the midpoint
    !! section is used for assumed constant axial and shear strain.
    class(beam_section), intent(in) :: section
    real(real64), intent(in) :: x(3,2), y(3), stiffness(3), shear(2)
        !! Reference geometry, local-y hint, [GJ, EI_y, EI_z], and [GA_y, GA_z].
    type(hyperdual), intent(in) :: u(:)
        !! Spatial nodal DOFs.
    type(hyperdual) :: w
        !! Total axial, shear, bending and torsional strain energy.
    real(real64) :: l0, frame(3,3)
    type(hyperdual) :: chord(3), q(4,2), mid(4), relative(3), cross(3), gamma(3)
    type(hyperdual) :: cos_half, norm, ratio2, k, factor
    integer :: i, j

    if (.not. all(stiffness >= 0.0d0) .or. .not. all(shear >= 0.0d0)) &
        error stop "Flexible beam stiffness cannot be negative"
    call beam_kinematics(section, x, y, u, l0, frame, chord, q)

    cos_half = hd_dot(q(:,1), q(:,2))
    if (.not. (abs(cos_half%x) > frame_tol)) error stop "Timoshenko relative rotation is at the pi branch cut"
    if (cos_half%x < 0.0d0) then
        cos_half = -cos_half
        do i = 1, 4
            q(i,2) = -q(i,2)
        end do
    end if
    norm = sqrt(2.0d0 + 2.0d0*cos_half)
    do i = 1, 4
        mid(i) = (q(i,1) + q(i,2)) / norm
    end do

    ! Vector part of conjugate(q1)*q2, expressed in the global reference basis.
    cross = hd_cross(q(2:4,1), q(2:4,2))
    do i = 1, 3
        relative(i) = q(1,1)*q(i+1,2) - q(1,2)*q(i+1,1) - cross(i)
    end do
    norm = hd_dot(relative, relative)
    if (norm%x < 1.0d-8) then
        ratio2 = norm / (cos_half*cos_half)
        factor = (2.0d0 / l0) * (1.0d0 - ratio2/3.0d0 + ratio2*ratio2/5.0d0) / cos_half
    else
        norm = sqrt(norm)
        factor = (2.0d0 / l0) * atan(norm / cos_half) / norm
    end if
    do i = 1, 3
        gamma(i) = hd_dot(chord, rotate_axis(mid, frame(:,i))) / l0
    end do
    w = section%area * section%material%strain_energy_density(gamma(1) - 1.0d0)
    do i = 1, 3
        k = 0.0d0
        do j = 1, 3
            k = k + relative(j)*frame(j,i)
        end do
        k = factor*k
        w = w + 0.5d0 * stiffness(i)*k*k
    end do
    do i = 1, 2
        w = w + 0.5d0 * shear(i)*gamma(i+1)*gamma(i+1)
    end do
    w = l0*w
end function

function rotation_quaternion(r) result(q)
    !! Maps a global rotation vector to a scalar-first unit quaternion. In the
    !! small-angle branch the series in angle squared preserves the Hessian at zero.
    type(hyperdual), intent(in) :: r(3)
        !! Rotation vector in radians.
    type(hyperdual) :: q(4)
        !! Unit quaternion, with its scalar component first.
    type(hyperdual) :: a2, half, scale
    integer :: i

    a2 = hd_dot(r, r)
    if (a2%x < 1.0d-8) then
        q(1) = 1.0d0 - a2/8.0d0 + a2*a2/384.0d0
        scale = 0.5d0 - a2/48.0d0 + a2*a2/3840.0d0
    else
        half = 0.5d0*sqrt(a2)
        q(1) = cos(half)
        scale = 0.5d0*sin(half)/half
    end if
    do i = 1, 3
        q(i+1) = scale*r(i)
    end do
end function

function rotate_axis(q, a) result(b)
    !! Rotates the real reference director `a` by the unit quaternion `q`.
    type(hyperdual), intent(in) :: q(4)
        !! Scalar-first unit quaternion.
    real(real64), intent(in) :: a(3)
        !! Director in the reference configuration.
    type(hyperdual) :: b(3)
        !! Rotated director.
    type(hyperdual) :: ah(3), cross(3), double_cross(3)
    integer :: i

    do i = 1, 3
        ah(i) = a(i)
    end do
    cross = hd_cross(q(2:4), ah)
    double_cross = hd_cross(q(2:4), cross)
    do i = 1, 3
        b(i) = a(i) + 2.0d0*(q(1)*cross(i) + double_cross(i))
    end do
end function

function hd_dot(a, b) result(c)
    !! Returns the dot product, retaining derivatives of both operands.
    type(hyperdual), intent(in) :: a(:), b(size(a))
    type(hyperdual) :: c
    integer :: i

    c = 0.0d0
    do i = 1, size(a)
        c = c + a(i)*b(i)
    end do
end function

function hd_cross(a, b) result(c)
    !! Returns the right-handed 3D cross product, retaining both operands' derivatives.
    type(hyperdual), intent(in) :: a(3), b(3)
    type(hyperdual) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
end function

end module
