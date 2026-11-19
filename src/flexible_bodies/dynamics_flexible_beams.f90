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
    !! kinematics and replaceable section energy laws.
    use iso_fortran_env, only : real64, int32
    use ieee_arithmetic, only : ieee_is_finite
    use autodiff
    use dynamics_helper, only : cross_product
    use dynamics_flexible_types, only : finite_element
    use dynamics_flexible_materials, only : flexible_material
    implicit none
    private
    public :: flexible_euler_bernoulli_2d, flexible_euler_bernoulli_3d
    public :: flexible_timoshenko_2d, flexible_timoshenko_3d
    public :: beam_section_geometry, beam_section_strain
    public :: beam_section_energy_law, quadratic_beam_section_energy

    real(real64), parameter :: frame_tol = 100.0d0 * epsilon(1.0d0)
        !! Tolerance for degenerate section frames and the relative rotation branch cut.

    type :: beam_section_geometry
        !! Reference section geometry in the beam's local axes.
        real(real64) :: area
            !! Positive cross-sectional area \(A\), in units of length squared.
        real(real64) :: second_moment_of_area(2)
            !! Second moments \([I_y, I_z]\), in units of length to the fourth power.
        real(real64) :: torsion_constant = 0.0d0
            !! Saint-Venant torsion constant \(J\), in the same units as the second moments.
        real(real64) :: effective_shear_area(2) = 0.0d0
            !! Effective areas along local [y, z], including shear correction.
    end type

    type :: beam_section_strain
        !! Local strain measures at one section evaluation point. In planar
        !! beams only extension, z curvature and (for Timoshenko) y shear are active.
        type(hyperdual) :: extension
            !! Dimensionless axial stretch minus one (projected onto the director for Timoshenko).
        type(hyperdual) :: shear(2)
            !! Transverse shear along local [y, z]; zero for Euler-Bernoulli.
        type(hyperdual) :: curvature(3)
            !! Twist and bending rates about local [x, y, z] per reference length.
    end type

    type, abstract :: beam_section_energy_law
        !! Conservative constitutive response of a beam section. Implementations
        !! preserve hyper-dual derivatives so the element can obtain forces and a tangent.
    contains
        procedure(section_energy_function), deferred :: energy_density
            !! Computes stored energy per unit reference length.
    end type

    type, extends(beam_section_energy_law) :: quadratic_beam_section_energy
        !! Default uncoupled section law. The axial contribution uses the
        !! material's energy law and need not be quadratic.
    contains
        procedure :: energy_density => quadratic_energy_density
            !! Evaluates the default section response per unit reference length.
    end type

    type(quadratic_beam_section_energy), parameter :: default_section_law = quadratic_beam_section_energy()

    abstract interface
        function section_energy_function(this, geometry, material, strain) result(w)
            !! Evaluates the stored energy at one section integration point.
            import :: beam_section_energy_law, beam_section_geometry, flexible_material, beam_section_strain, hyperdual
            class(beam_section_energy_law), intent(in) :: this
                !! Section law and any parameters it carries.
            type(beam_section_geometry), intent(in) :: geometry
                !! Section geometry in the reference configuration.
            class(flexible_material), intent(in) :: material
                !! Constitutive material assigned to the element.
            type(beam_section_strain), intent(in) :: strain
                !! Local generalized strains, including hyper-dual derivative components.
            type(hyperdual) :: w
                !! Stored energy per unit reference length.
        end function
    end interface

    type, abstract, extends(finite_element) :: beam_section
        !! Shared section inputs for planar and spatial beam elements.
        real(real64) :: area = 1.0d0
            !! Reference cross-sectional area.
        class(flexible_material), allocatable :: material
            !! Constitutive material passed to the section law; must be allocated before evaluation.
        class(beam_section_energy_law), allocatable :: section_law
            !! If unallocated, uses quadratic_beam_section_energy.
    end type

    type, abstract, extends(beam_section) :: planar_beam
        !! Base for two-node beams confined to the global xy plane.
        real(real64) :: reference_coordinates(2,2) = 0.0d0
            !! One column per node. DOFs are `[ux, uy, theta_z]` per node,
            !! with positive angles measured counterclockwise in radians.
        real(real64) :: second_moment_of_area = 1.0d0
            !! Reference second moment of area I_z about the out-of-plane axis.
    contains
        procedure :: get_dof_count => planar_beam_dof_count
    end type

    type, abstract, extends(beam_section) :: spatial_beam
        !! Base for two-node beams with independent translations and finite rotations in 3D.
        real(real64) :: reference_coordinates(3,2) = 0.0d0
            !! One column per node. DOFs are `[ux, uy, uz, rx, ry, rz]` per node.
            !! The last three are a global rotation vector in radians, applied
            !! through the exponential map to the reference section frame.
        real(real64) :: reference_y(3) = [0.0d0, 1.0d0, 0.0d0]
            !! Direction projected onto the reference cross section to set its
            !! local y axis; it must not be parallel to the reference chord.
        real(real64) :: second_moment_of_area(2) = 1.0d0
            !! Reference second moments of area [I_y, I_z].
        real(real64) :: torsion_constant = 1.0d0
            !! Reference Saint-Venant torsion constant J.
    contains
        procedure :: get_dof_count => spatial_beam_dof_count
    end type

    type, extends(planar_beam) :: flexible_euler_bernoulli_2d
        !! Shear-free planar beam with a cubic centerline.
        !! Section energy is integrated with four Gauss points.
    contains
        procedure :: strain_energy => euler_bernoulli_2d_energy
    end type

    type, extends(spatial_beam) :: flexible_euler_bernoulli_3d
        !! Shear-free spatial beam with a cubic centerline and section twist.
        !! The projected section frame must remain nonsingular at all four Gauss points.
    contains
        procedure :: strain_energy => euler_bernoulli_3d_energy
        procedure :: mass_matrix => euler_bernoulli_3d_mass
    end type

    type, extends(planar_beam) :: flexible_timoshenko_2d
        !! Extensible planar beam with independent section rotations and shear.
        !! Extension and shear are evaluated at the midpoint.
        real(real64) :: effective_shear_area = 1.0d0
            !! Reference shear area along local y, including the shear correction factor.
    contains
        procedure :: strain_energy => timoshenko_2d_energy
    end type

    type, extends(spatial_beam) :: flexible_timoshenko_3d
        !! Extensible spatial beam with independent section rotations, shear and twist.
        !! Shortest-arc rotation interpolation is undefined at a relative rotation of pi.
        real(real64) :: effective_shear_area(2) = 1.0d0
            !! Reference shear areas along local [y, z], including shear correction.
    contains
        procedure :: strain_energy => timoshenko_3d_energy
        procedure :: mass_matrix => timoshenko_3d_mass
    end type

contains

function euler_bernoulli_3d_mass(this) result(m)
    !! Reference consistent mass including section rotary inertia.
    class(flexible_euler_bernoulli_3d), intent(in) :: this
    real(real64), allocatable :: m(:,:)

    m = spatial_beam_mass(this, .true.)
end function

function timoshenko_3d_mass(this) result(m)
    !! Reference consistent mass with independent linear rotation interpolation.
    class(flexible_timoshenko_3d), intent(in) :: this
    real(real64), allocatable :: m(:,:)

    m = spatial_beam_mass(this, .false.)
end function

function spatial_beam_mass(this, shear_free) result(m)
    !! Integrates the reference velocity and angular velocity shape functions.
    !! In the shear-free case, bending rotations are slopes of the Hermite
    !! centerline. The rotary inertia about the chord uses Iy+Iz, not J.
    class(spatial_beam), intent(in) :: this
    logical, intent(in) :: shear_free
    real(real64) :: m(12,12), frame(3,3), length, x, n(2), h(4), dh(4)
    real(real64) :: b(6,12), weighted(6,12), inertia(6)
    type(hyperdual) :: u(12), chord(3), q(4,2)
    real(real64), parameter :: points(4) = [0.0694318442029737d0, 0.330009478207572d0, &
        0.669990521792428d0, 0.930568155797026d0]
    real(real64), parameter :: weights(4) = [0.173927422568727d0, 0.326072577431273d0, &
        0.326072577431273d0, 0.173927422568727d0]
    integer :: i, j, k

    call validate_section(this, beam_section_geometry(this%area, this%second_moment_of_area, this%torsion_constant))
    if (.not. ieee_is_finite(this%material%density) .or. this%material%density < 0.0d0) &
        error stop "Invalid flexible beam density"
    do i = 1, 12
        u(i) = 0.0d0
    end do
    call beam_kinematics(this%reference_coordinates, this%reference_y, u, length, frame, chord, q)
    inertia = [this%area, this%area, this%area, sum(this%second_moment_of_area), this%second_moment_of_area]
    m = 0.0d0
    do i = 1, 4
        x = points(i)
        n = [1.0d0-x, x]
        b = 0.0d0
        do j = 1, 6
            b(j,j) = n(1)
            b(j,j+6) = n(2)
        end do
        if (shear_free) then
            h = [1.0d0-3*x**2+2*x**3, length*x*(1-x)**2, 3*x**2-2*x**3, length*x**2*(x-1)]
            dh = [6*x*(x-1)/length, 1-4*x+3*x**2, 6*x*(1-x)/length, 3*x**2-2*x]
            b(2,[2,6,8,12]) = h
            b(3,[3,5,9,11]) = h * [1,-1,1,-1]
            b(5,[3,5,9,11]) = dh * [-1,1,-1,1]
            b(6,[2,6,8,12]) = dh
        end if
        do k = 1, 12, 3
            b(:,k:k+2) = matmul(b(:,k:k+2), transpose(frame))
        end do
        weighted = b * spread(inertia, 2, 12)
        m = m + weights(i) * matmul(transpose(b), weighted)
    end do
    m = m * length * this%material%density
end function

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
        !! Planar beam providing reference coordinates, section geometry, and its energy law.
    type(hyperdual), intent(in) :: u(:)
        !! Six nodal translations and angles in node order.
    type(hyperdual) :: w
        !! Total strain energy.
    type(hyperdual) :: u3(12)
    real(real64) :: x(3,2), y(3)

    call embed_planar(this%reference_coordinates, u, x, y, u3)
    w = euler_bernoulli_energy(this, x, y, &
        beam_section_geometry(this%area, [0.0d0, this%second_moment_of_area]), u3)
end function

function euler_bernoulli_3d_energy(this, u) result(w)
    !! Evaluates the spatial axial, bending and torsional energy.
    class(flexible_euler_bernoulli_3d), intent(in) :: this
        !! Spatial beam providing reference geometry and its section energy law.
    type(hyperdual), intent(in) :: u(:)
        !! Twelve nodal translations and rotation-vector components in node order.
    type(hyperdual) :: w
        !! Total strain energy.

    w = euler_bernoulli_energy(this, this%reference_coordinates, this%reference_y, &
        beam_section_geometry(this%area, this%second_moment_of_area, this%torsion_constant), u)
end function

function timoshenko_2d_energy(this, u) result(w)
    !! Evaluates the planar axial, bending and shear energy.
    class(flexible_timoshenko_2d), intent(in) :: this
        !! Planar beam providing reference coordinates, section geometry, and its energy law.
    type(hyperdual), intent(in) :: u(:)
        !! Six nodal translations and angles in node order.
    type(hyperdual) :: w
        !! Total strain energy.
    type(hyperdual) :: u3(12)
    real(real64) :: x(3,2), y(3)

    call embed_planar(this%reference_coordinates, u, x, y, u3)
    w = timoshenko_energy(this, x, y, beam_section_geometry(this%area, [0.0d0, this%second_moment_of_area], &
        effective_shear_area = [this%effective_shear_area, 0.0d0]), u3)
end function

function timoshenko_3d_energy(this, u) result(w)
    !! Evaluates the spatial axial, bending, shear and torsional energy.
    class(flexible_timoshenko_3d), intent(in) :: this
        !! Spatial beam providing reference geometry and its section energy law.
    type(hyperdual), intent(in) :: u(:)
        !! Twelve nodal translations and rotation-vector components in node order.
    type(hyperdual) :: w
        !! Total strain energy.

    w = timoshenko_energy(this, this%reference_coordinates, this%reference_y, &
        beam_section_geometry(this%area, this%second_moment_of_area, this%torsion_constant, this%effective_shear_area), u)
end function

subroutine embed_planar(x, u, x3, y, u3)
    !! Embeds the two planar reference nodes and six active DOFs in the xy
    !! plane of the spatial formulation. The reference local y axis is left
    !! of the chord so local z always coincides with global z.
    !! Stops with an error unless exactly six planar DOFs are supplied.
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

subroutine beam_kinematics(x, y, u, l0, frame, chord, q)
    !! Validates the reference frame, then returns the reference
    !! length/frame, current chord, and the two global rotation quaternions.
    !! Requires twelve DOFs, a nonzero reference chord, and a nonzero local-y
    !! hint that is not parallel to the chord; violations stop with an error.
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

subroutine validate_section(section, geometry)
    !! Stops with an error if the material is unallocated, the area is not
    !! finite and positive, or any other section property is not finite and nonnegative.
    class(beam_section), intent(in) :: section
        !! Beam owning the constitutive material.
    type(beam_section_geometry), intent(in) :: geometry
        !! Reference geometry supplied to the section law.
    real(real64) :: moments(3), shear_area(2)

    if (.not. allocated(section%material)) error stop "Flexible beam requires a material"
    if (.not. ieee_is_finite(geometry%area) .or. .not. (geometry%area > 0.0d0)) &
        error stop "Flexible beam requires finite positive area"
    moments = [geometry%torsion_constant, geometry%second_moment_of_area]
    if (.not. all(ieee_is_finite(moments)) .or. .not. all(moments >= 0.0d0)) &
        error stop "Flexible beam requires finite nonnegative second moments and torsion constant"
    shear_area = geometry%effective_shear_area
    if (.not. all(ieee_is_finite(shear_area)) .or. .not. all(shear_area >= 0.0d0)) &
        error stop "Flexible beam requires finite nonnegative effective shear areas"
end subroutine

function section_energy_density(section, geometry, strain) result(w)
    !! Dispatches to the assigned section law or the built-in default.
    class(beam_section), intent(in) :: section
        !! Beam with an allocated material and an optional section law.
    type(beam_section_geometry), intent(in) :: geometry
        !! Validated reference section geometry.
    type(beam_section_strain), intent(in) :: strain
        !! Hyper-dual strains in the local section axes.
    type(hyperdual) :: w
        !! Stored energy per unit reference length.

    if (allocated(section%section_law)) then
        w = section%section_law%energy_density(geometry, section%material, strain)
    else
        w = default_section_law%energy_density(geometry, section%material, strain)
    end if
end function

function quadratic_energy_density(this, geometry, material, strain) result(w)
    !! Multiplies the material's axial energy density by the reference area
    !! and adds independent quadratic terms weighted by \(GJ\), \(EI_y\),
    !! \(EI_z\), and the effective shear stiffnesses \(GA_s\).
    !! Stops with an error for nonpositive or nonfinite moduli, or nonfinite stiffnesses.
    class(quadratic_beam_section_energy), intent(in) :: this
        !! Stateless default section law.
    type(beam_section_geometry), intent(in) :: geometry
        !! Reference area, second moments, torsion constant, and effective shear areas.
    class(flexible_material), intent(in) :: material
        !! Material providing the axial energy law and constant elastic moduli.
    type(beam_section_strain), intent(in) :: strain
        !! Hyper-dual axial, shear, and curvature measures.
    type(hyperdual) :: w
        !! Stored energy per unit reference length.
    real(real64) :: e, g, stiffness(3), shear(2)
    integer :: i

    call material%get_elastic_moduli(e, g)
    if (.not. all(ieee_is_finite([e, g])) .or. .not. all([e, g] > 0.0d0)) &
        error stop "Flexible beam requires finite positive Young's and shear moduli"
    stiffness = [g, e, e] * [geometry%torsion_constant, geometry%second_moment_of_area]
    shear = g * geometry%effective_shear_area
    if (.not. all(ieee_is_finite(stiffness))) error stop "Flexible beam curvature stiffness must be finite"
    if (.not. all(ieee_is_finite(shear))) error stop "Flexible beam shear stiffness must be finite"

    w = geometry%area * material%strain_energy_density(strain%extension)
    do i = 1, 3
        w = w + 0.5d0 * stiffness(i) * strain%curvature(i) * strain%curvature(i)
    end do
    do i = 1, 2
        w = w + 0.5d0 * shear(i) * strain%shear(i) * strain%shear(i)
    end do
end function

function euler_bernoulli_energy(section, x, y, geometry, u) result(w)
    !! Integrates the shear-free energy with four Gauss points on the reference
    !! interval. The nodal section axes set the Hermite tangent directions;
    !! the current chord length scales both end derivatives. The local y axis
    !! is interpolated then projected normal to the current centerline.
    !! Stops with an error if the current chord vanishes or if the tangent
    !! or projected section frame is singular at a quadrature point.
    class(beam_section), intent(in) :: section
        !! Beam providing the material and section law.
    real(real64), intent(in) :: x(3,2), y(3)
        !! Reference nodal coordinates (columns) and local-y orientation hint.
    type(beam_section_geometry), intent(in) :: geometry
        !! Reference section geometry; effective shear areas are zero.
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
    type(hyperdual) :: dr(3), ddr(3), t(3), dt(3), guide(3), dguide(3), d2(3), d3(3)
    type(beam_section_strain) :: strain
    integer :: i, j

    call validate_section(section, geometry)
    call beam_kinematics(x, y, u, l0, frame, chord, q)
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

    strain%shear(1) = 0.0d0
    strain%shear(2) = 0.0d0
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
        strain%curvature(1) = (hd_dot(d3, dguide) - projection*hd_dot(d3, dt)) / width
        strain%curvature(2) = -hd_dot(d3, dt)
        strain%curvature(3) = hd_dot(d2, dt)
        strain%extension = speed - 1.0d0
        w = w + weight(j) * l0 * section_energy_density(section, geometry, strain)
    end do
end function

function timoshenko_energy(section, x, y, geometry, u) result(w)
    !! Uses a linear centerline and shortest-arc section interpolation. The
    !! relative rotation gives constant material curvature; the midpoint
    !! section is used for assumed constant axial and shear strain.
    !! Stops with an error at the relative rotation branch cut near pi.
    class(beam_section), intent(in) :: section
        !! Beam providing the material and section law.
    real(real64), intent(in) :: x(3,2), y(3)
        !! Reference nodal coordinates (columns) and local-y orientation hint.
    type(beam_section_geometry), intent(in) :: geometry
        !! Reference section geometry, including effective shear areas.
    type(hyperdual), intent(in) :: u(:)
        !! Spatial nodal DOFs.
    type(hyperdual) :: w
        !! Total axial, shear, bending and torsional strain energy.
    real(real64) :: l0, frame(3,3)
    type(hyperdual) :: chord(3), q(4,2), mid(4), relative(3), cross(3)
    type(hyperdual) :: cos_half, norm, ratio2, k, factor
    type(beam_section_strain) :: strain
    integer :: i, j

    call validate_section(section, geometry)
    call beam_kinematics(x, y, u, l0, frame, chord, q)

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
    strain%extension = hd_dot(chord, rotate_axis(mid, frame(:,1))) / l0 - 1.0d0
    do i = 1, 2
        strain%shear(i) = hd_dot(chord, rotate_axis(mid, frame(:,i+1))) / l0
    end do
    do i = 1, 3
        k = 0.0d0
        do j = 1, 3
            k = k + relative(j)*frame(j,i)
        end do
        strain%curvature(i) = factor*k
    end do
    w = l0 * section_energy_density(section, geometry, strain)
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
        !! Equal-length vector operands.
    type(hyperdual) :: c
        !! Scalar dot product.
    integer :: i

    c = 0.0d0
    do i = 1, size(a)
        c = c + a(i)*b(i)
    end do
end function

function hd_cross(a, b) result(c)
    !! Returns the right-handed 3D cross product, retaining both operands' derivatives.
    type(hyperdual), intent(in) :: a(3), b(3)
        !! Vector operands in a common right-handed basis.
    type(hyperdual) :: c(3)
        !! Cross product \(a \times b\) in the same basis.

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
end function

end module
