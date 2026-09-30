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

! References:
! - J.-L. Batoz, K.-J. Bathe, L.-W. Ho, "A study of three-node triangular
!   plate bending elements," Int. J. Numer. Meth. Engng., 15, 1771-1812, 1980.
! - E. N. Dvorkin, K.-J. Bathe, "A continuum mechanics based four-node shell
!   element for general nonlinear analysis," Eng. Comput., 1, 77-88, 1984.
! - T. J. R. Hughes, F. Brezzi, "On drilling degrees of freedom," Comput.
!   Methods Appl. Mech. Engrg., 72, 105-121, 1989.
! - R. D. Cook, D. S. Malkus, M. E. Plesha, R. J. Witt, "Concepts and
!   Applications of Finite Element Analysis," 4th ed., Wiley, 2002.

module dynamics_shell_elements
    !! Flat, three-dimensional shell elements with six degrees of freedom per
    !! node.
    !!
    !! Each element is formulated in a local, planar coordinate system by
    !! superposing a membrane (plane-stress) element, a plate-bending element,
    !! and a drilling-rotation stabilization term.  The resulting local
    !! matrices are then transformed into the global coordinate system.
    !!
    !! The nodal degrees of freedom are ordered as
    !! $$ \boldsymbol{u}_i = \begin{bmatrix} u_i & v_i & w_i & \theta_{x,i} &
    !! \theta_{y,i} & \theta_{z,i} \end{bmatrix}^T, $$
    !! where the rotations follow the right-hand rule about each axis.  The
    !! element displacement vector stacks the nodal vectors in local node
    !! order.
    !!
    !! **Local Coordinate System**
    !!
    !! The local \(z\)-axis is the unit normal of the element plane computed
    !! from the node ordering (right-hand rule).  The local \(x\)-axis is the
    !! projection of the vector from node 1 to node 2 onto the element plane,
    !! and the local \(y\)-axis completes the right-handed system.  Nodes are
    !! assumed to lie in a common plane; warped quadrilaterals are projected
    !! onto their mean plane.
    !!
    !! **Kinematics**
    !!
    !! Within the local system, the displacement through the thickness is
    !! \(u = u_0 + z\theta_y\) and \(v = v_0 - z\theta_x\).  The generalized
    !! strain vector is
    !! $$ \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} &
    !! \varepsilon_{yy} & \gamma_{xy} & \kappa_{xx} & \kappa_{yy} &
    !! \kappa_{xy} & \gamma_{xz} & \gamma_{yz} \end{bmatrix}^T, $$
    !! with
    !! $$ \kappa_{xx} = \frac{\partial \theta_y}{\partial x}, \quad
    !! \kappa_{yy} = -\frac{\partial \theta_x}{\partial y}, \quad
    !! \kappa_{xy} = \frac{\partial \theta_y}{\partial y} -
    !! \frac{\partial \theta_x}{\partial x}, $$
    !! $$ \gamma_{xz} = \frac{\partial w}{\partial x} + \theta_y, \quad
    !! \gamma_{yz} = \frac{\partial w}{\partial y} - \theta_x. $$
    !!
    !! The corresponding stress-resultant vector is
    !! $$ \boldsymbol{\sigma} = \begin{bmatrix} N_{xx} & N_{yy} & N_{xy} &
    !! M_{xx} & M_{yy} & M_{xy} & Q_x & Q_y \end{bmatrix}^T = D
    !! \boldsymbol{\varepsilon}, $$
    !! where \(N\) are membrane forces per unit length, \(M\) are moments per
    !! unit length, and \(Q\) are transverse shear forces per unit length.
    !!
    !! **Drilling Rotation**
    !!
    !! The in-plane rotation \(\theta_z\) is not part of classical shell
    !! theory.  It is stabilized with the Hughes-Brezzi penalty
    !! $$ \Pi_d = \frac{1}{2}\int_A \alpha G t \left( \theta_z -
    !! \frac{1}{2}\left( \frac{\partial v}{\partial x} -
    !! \frac{\partial u}{\partial y} \right) \right)^2 dA, $$
    !! which vanishes for rigid-body motion and prevents a singular stiffness
    !! matrix for co-planar meshes.  The factor \(\alpha\) is the element's
    !! drilling_factor.
    !!
    !! **Mass Matrix**
    !!
    !! The consistent mass matrix is computed from the linear (triangle) or
    !! bilinear (quadrilateral) interpolation of all six nodal quantities with
    !! translational inertia \(\rho t\) and rotary inertia
    !! \(\rho t^3 / 12\).  The rotary inertia is also assigned to the drilling
    !! rotation to keep the mass matrix positive definite.
    use iso_fortran_env, only : int32, real64
    use dynamics_error_handling
    use dynamics_helper, only : cross_product
    use dynamics_structural, only : node, material, element
    implicit none
    private
    public :: shell_element
    public :: triangular_shell_element
    public :: rectangular_shell_element

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type, extends(element), abstract :: shell_element
        !! Defines a flat shell element with six degrees of freedom per node
        !! (three translations and three rotations).
        real(real64) :: thickness
            !! The shell thickness.
        real(real64) :: drilling_factor = 1.0d-3
            !! The nondimensional drilling-rotation penalty factor
            !! \(\alpha\).  The drilling stiffness is \(\alpha G t\), where
            !! \(G\) is the shear modulus.
        real(real64) :: shear_correction = 5.0d0 / 6.0d0
            !! The transverse shear correction factor.  This value is only
            !! used by elements that include transverse shear deformation.
    contains
        procedure(shell_natural_gradient), deferred, public, pass :: &
            shape_function_natural_gradient
        procedure(shell_integration_rule), deferred, public, pass :: &
            integration_rule
        procedure, public :: get_dimensionality => shl_dimensionality
        procedure, public :: get_dof_per_node => shl_dof_per_node
        procedure, public :: constitutive_matrix => shl_constitutive_matrix
        procedure, public :: jacobian => shl_jacobian
        procedure, public :: shape_function_matrix => shl_shape_function_matrix
        procedure, public :: shape_function_gradient => &
            shl_shape_function_gradient
        procedure, public :: local_frame => shl_local_frame
        procedure, public :: local_coordinates => shl_local_coordinates
        procedure, public :: rotation_matrix => shl_rotation_matrix
        procedure, public :: area => shl_area
        procedure, public :: stiffness_matrix => shl_stiffness_matrix
        procedure, public :: mass_matrix => shl_mass_matrix
        procedure, public :: external_force_vector => shl_ext_force_vector
        procedure, public :: strain => shl_strain
        procedure, public :: stress => shl_stress
    end type

! ------------------------------------------------------------------------------
    type, extends(shell_element) :: triangular_shell_element
        !! Defines a three-node flat shell element.
        !!
        !! The membrane behavior is modeled by the constant-strain triangle
        !! (CST) and the bending behavior by the Discrete Kirchhoff Triangle
        !! (DKT) of Batoz, Bathe, and Ho.  The DKT element enforces the
        !! Kirchhoff hypothesis at discrete points; therefore, transverse
        !! shear deformation is neglected and the transverse shear strain and
        !! force components are always zero.  This element is appropriate
        !! for thin shells.
        !!
        !! The natural coordinates \((\xi, \eta)\) are the area coordinates
        !! of nodes 2 and 3, so nodes 1, 2, and 3 are located at \((0,0)\),
        !! \((1,0)\), and \((0,1)\), respectively.  Nodes should be ordered
        !! counter-clockwise when viewed from the positive local \(z\)-axis;
        !! the local \(z\)-axis is defined by this ordering.
        type(node), dimension(3) :: nodes
            !! The element nodes.
    contains
        procedure, public :: get_node_count => tri_get_node_count
        procedure, public :: get_node => tri_get_node
        procedure, public :: get_node_natural_coordinates => &
            tri_get_node_natural_coordinates
        procedure, public :: evaluate_shape_function => tri_shape_function
        procedure, public :: shape_function_natural_gradient => &
            tri_shape_function_natural_gradient
        procedure, public :: integration_rule => tri_integration_rule
        procedure, public :: strain_displacement_matrix => &
            tri_strain_disp_matrix
    end type

    interface triangular_shell_element
        module procedure :: tri_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(shell_element) :: rectangular_shell_element
        !! Defines a four-node flat shell element.
        !!
        !! The membrane behavior is modeled by the bilinear isoparametric
        !! quadrilateral and the bending behavior by the Mindlin-Reissner
        !! MITC4 formulation of Bathe and Dvorkin, which uses assumed
        !! transverse shear strains tied at the edge midpoints to eliminate
        !! shear locking.  Both thick and thin shells are supported.
        !! Although nominally rectangular, the isoparametric formulation is
        !! valid for any flat, convex quadrilateral.
        !!
        !! The natural coordinates \((r, s)\) span \([-1, 1]^2\), with nodes
        !! 1 through 4 located at \((-1,-1)\), \((1,-1)\), \((1,1)\), and
        !! \((-1,1)\), respectively.  Nodes should be ordered
        !! counter-clockwise when viewed from the positive local \(z\)-axis;
        !! the local \(z\)-axis is defined by this ordering.
        type(node), dimension(4) :: nodes
            !! The element nodes.
    contains
        procedure, public :: get_node_count => quad_get_node_count
        procedure, public :: get_node => quad_get_node
        procedure, public :: get_node_natural_coordinates => &
            quad_get_node_natural_coordinates
        procedure, public :: evaluate_shape_function => quad_shape_function
        procedure, public :: shape_function_natural_gradient => &
            quad_shape_function_natural_gradient
        procedure, public :: integration_rule => quad_integration_rule
        procedure, public :: strain_displacement_matrix => &
            quad_strain_disp_matrix
    end type

    interface rectangular_shell_element
        module procedure :: quad_init
    end interface

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        pure function shell_natural_gradient(this, s) result(rst)
            !! Defines the signature of a routine returning the derivatives of
            !! the element shape functions with respect to the natural
            !! coordinates.
            import :: shell_element, real64
            class(shell_element), intent(in) :: this
                !! The shell_element object.
            real(real64), intent(in), dimension(:) :: s
                !! The natural coordinates at which to evaluate the
                !! derivatives.
            real(real64), allocatable, dimension(:,:) :: rst
                !! A 2-by-N matrix, where N is the number of nodes, whose
                !! first row contains the derivatives with respect to the
                !! first natural coordinate and whose second row contains the
                !! derivatives with respect to the second natural coordinate.
        end function

        pure subroutine shell_integration_rule(this, pts, wts)
            !! Defines the signature of a routine returning the numerical
            !! integration rule used by the element.  An integral over the
            !! element is evaluated as
            !! \(\int_A f\,dA \approx \sum_i w_i f(\boldsymbol{s}_i)
            !! \det J(\boldsymbol{s}_i)\).
            import :: shell_element, real64
            class(shell_element), intent(in) :: this
                !! The shell_element object.
            real(real64), allocatable, intent(out), dimension(:,:) :: pts
                !! A 2-by-M matrix containing the natural coordinates of the
                !! M integration points.
            real(real64), allocatable, intent(out), dimension(:) :: wts
                !! The M integration weights.
        end subroutine
    end interface

contains
! ******************************************************************************
! PRIVATE HELPERS
! ------------------------------------------------------------------------------
pure function det2(x) result(rst)
    ! Computes the determinant of a 2-by-2 matrix.
    real(real64), intent(in), dimension(2,2) :: x
    real(real64) :: rst
    rst = x(1,1) * x(2,2) - x(1,2) * x(2,1)
end function

! ------------------------------------------------------------------------------
pure function node_positions(elem) result(rst)
    ! Returns a 3-by-N matrix of the global nodal positions.
    class(shell_element), intent(in) :: elem
    real(real64), allocatable, dimension(:,:) :: rst

    integer(int32) :: i, n
    type(node) :: nd

    n = elem%get_node_count()
    allocate(rst(3, n))
    do i = 1, n
        nd = elem%get_node(i)
        rst(:,i) = [nd%x, nd%y, nd%z]
    end do
end function

! ------------------------------------------------------------------------------
pure function membrane_strain_rows(dndx) result(rst)
    ! Builds the 3-by-6N membrane strain-displacement rows from the 2-by-N
    ! shape function gradient in local coordinates.
    real(real64), intent(in), dimension(:,:) :: dndx
    real(real64), allocatable, dimension(:,:) :: rst

    integer(int32) :: i, c

    allocate(rst(3, 6 * size(dndx, 2)), source = 0.0d0)
    do i = 1, size(dndx, 2)
        c = 6 * (i - 1)
        rst(1,c+1) = dndx(1,i)
        rst(2,c+2) = dndx(2,i)
        rst(3,c+1) = dndx(2,i)
        rst(3,c+2) = dndx(1,i)
    end do
end function

! ------------------------------------------------------------------------------
pure function drilling_row(elem, s) result(rst)
    ! Builds the 6N-element operator giving theta_z - 0.5 * (dv/dx - du/dy)
    ! at the natural coordinate s.
    class(shell_element), intent(in) :: elem
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:) :: rst

    integer(int32) :: i, c, n
    real(real64), allocatable, dimension(:,:) :: dndx

    n = elem%get_node_count()
    dndx = elem%shape_function_gradient(s)
    allocate(rst(6 * n), source = 0.0d0)
    do i = 1, n
        c = 6 * (i - 1)
        rst(c+1) = 0.5d0 * dndx(2,i)
        rst(c+2) = -0.5d0 * dndx(1,i)
        rst(c+6) = elem%evaluate_shape_function(i, s)
    end do
end function

! ******************************************************************************
! SHELL_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function shl_dimensionality(this) result(rst)
    !! Gets the dimensionality of the element.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    integer(int32) :: rst
        !! The dimensionality (always 3).
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function shl_dof_per_node(this) result(rst)
    !! Gets the number of degrees of freedom per node.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    integer(int32) :: rst
        !! The number of degrees of freedom per node (always 6).
    rst = 6
end function

! ------------------------------------------------------------------------------
pure function shl_constitutive_matrix(this) result(rst)
    !! Computes the 8-by-8 constitutive matrix relating the generalized
    !! strains to the stress resultants.
    !! $$ D = \begin{bmatrix} t C & 0 & 0 \\ 0 & \frac{t^3}{12} C & 0 \\
    !! 0 & 0 & k G t I \end{bmatrix}, \quad
    !! C = \frac{E}{1 - \nu^2} \begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\
    !! 0 & 0 & \frac{1 - \nu}{2} \end{bmatrix}, $$
    !! where \(k\) is the shear correction factor and
    !! \(G = E / (2 (1 + \nu))\).
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 8-by-8 constitutive matrix.

    ! Local Variables
    real(real64) :: e, nu, t, g, c(3,3)

    ! Initialization
    e = this%material%modulus
    nu = this%material%poissons_ratio
    t = this%thickness
    g = e / (2.0d0 * (1.0d0 + nu))

    ! Plane-stress material matrix
    c = reshape([1.0d0, nu, 0.0d0, nu, 1.0d0, 0.0d0, &
        0.0d0, 0.0d0, 0.5d0 * (1.0d0 - nu)], [3, 3])
    c = (e / (1.0d0 - nu**2)) * c

    ! Assemble the membrane, bending, and transverse shear blocks
    allocate(rst(8, 8), source = 0.0d0)
    rst(1:3,1:3) = t * c
    rst(4:6,4:6) = (t**3 / 12.0d0) * c
    rst(7,7) = this%shear_correction * g * t
    rst(8,8) = rst(7,7)
end function

! ------------------------------------------------------------------------------
pure function shl_local_frame(this) result(rst)
    !! Computes the direction cosine matrix of the element's local
    !! coordinate system.  The rows of the matrix are the local \(x\), \(y\),
    !! and \(z\) unit vectors expressed in global coordinates; therefore,
    !! \(\boldsymbol{x}_{local} = R \boldsymbol{x}_{global}\).
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 3-by-3 direction cosine matrix.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: nrm, scale, e1(3), e2(3), e3(3), d(3)
    real(real64), allocatable, dimension(:,:) :: p

    ! Nodal positions relative to node 1
    p = node_positions(this)
    n = size(p, 2)
    do i = n, 1, -1
        p(:,i) = p(:,i) - p(:,1)
    end do
    scale = maxval(norm2(p, dim = 1))
    if (scale <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR

    ! Newell's method yields the polygon normal (exact for planar polygons
    ! and a best-fit normal for mildly warped quadrilaterals).
    e3 = 0.0d0
    do i = 2, n - 1
        e3 = e3 + cross_product(p(:,i), p(:,i+1))
    end do
    nrm = norm2(e3)
    if (nrm <= sqrt(epsilon(nrm)) * scale**2) &
        error stop DYN_INVALID_INPUT_ERROR
    e3 = e3 / nrm

    ! Local x-axis: the node 1 to node 2 vector projected onto the plane
    d = p(:,2) - dot_product(p(:,2), e3) * e3
    nrm = norm2(d)
    if (nrm <= sqrt(epsilon(nrm)) * scale) error stop DYN_INVALID_INPUT_ERROR
    e1 = d / nrm
    e2 = cross_product(e3, e1)

    allocate(rst(3, 3))
    rst(1,:) = e1
    rst(2,:) = e2
    rst(3,:) = e3
end function

! ------------------------------------------------------------------------------
pure function shl_local_coordinates(this) result(rst)
    !! Computes the in-plane nodal coordinates in the element's local
    !! coordinate system.  Node 1 is located at the local origin.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An N-by-2 matrix, where N is the number of nodes, containing the
        !! local \(x\) and \(y\) coordinates of each node.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: q(3)
    real(real64), allocatable, dimension(:,:) :: r, p

    ! Process
    r = this%local_frame()
    p = node_positions(this)
    n = size(p, 2)
    allocate(rst(n, 2))
    do i = 1, n
        q = p(:,i) - p(:,1)
        rst(i,1) = dot_product(r(1,:), q)
        rst(i,2) = dot_product(r(2,:), q)
    end do
end function

! ------------------------------------------------------------------------------
pure function shl_rotation_matrix(this) result(rst)
    !! Computes the transformation matrix \(T\) relating the global element
    !! displacement vector to the local element displacement vector such
    !! that \(\boldsymbol{u}_{local} = T \boldsymbol{u}_{global}\).  The
    !! matrix is block diagonal with the 3-by-3 direction cosine matrix
    !! repeated for the translations and rotations of each node.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 6N-by-6N transformation matrix, where N is the number of
        !! nodes.

    ! Local Variables
    integer(int32) :: i, n
    real(real64), allocatable, dimension(:,:) :: r

    ! Process
    r = this%local_frame()
    n = 2 * this%get_node_count()
    allocate(rst(3 * n, 3 * n), source = 0.0d0)
    do i = 1, n
        rst(3*i-2:3*i,3*i-2:3*i) = r
    end do
end function

! ------------------------------------------------------------------------------
pure function shl_jacobian(this, s) result(rst)
    !! Computes the 2-by-2 Jacobian matrix of the mapping from natural to
    !! local coordinates.
    !! $$ J = \begin{bmatrix} \frac{\partial x}{\partial s_1} &
    !! \frac{\partial y}{\partial s_1} \\ \frac{\partial x}{\partial s_2} &
    !! \frac{\partial y}{\partial s_2} \end{bmatrix} $$
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates at which to evaluate the Jacobian.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 2-by-2 Jacobian matrix.
    rst = matmul(this%shape_function_natural_gradient(s), &
        this%local_coordinates())
end function

! ------------------------------------------------------------------------------
pure function shl_shape_function_gradient(this, s) result(rst)
    !! Computes the derivatives of the shape functions with respect to the
    !! local \(x\) and \(y\) coordinates,
    !! \(\nabla_{xy} N = J^{-1} \nabla_{s} N\).
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates at which to evaluate the derivatives.
    real(real64), allocatable, dimension(:,:) :: rst
        !! A 2-by-N matrix, where N is the number of nodes, containing
        !! \(\partial N_i / \partial x\) in the first row and
        !! \(\partial N_i / \partial y\) in the second row.

    ! Local Variables
    real(real64) :: jac(2,2), jinv(2,2), dj

    ! Process
    jac = this%jacobian(s)
    dj = det2(jac)
    ! A non-positive determinant indicates an inverted or non-convex element
    if (dj <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    jinv = reshape([jac(2,2), -jac(2,1), -jac(1,2), jac(1,1)], [2, 2]) / dj
    rst = matmul(jinv, this%shape_function_natural_gradient(s))
end function

! ------------------------------------------------------------------------------
pure function shl_shape_function_matrix(this, s) result(rst)
    !! Computes the 6-by-6N shape function matrix interpolating the local
    !! nodal quantities \([u, v, w, \theta_x, \theta_y, \theta_z]\) with the
    !! element's linear (triangle) or bilinear (quadrilateral) shape
    !! functions.  This interpolation is used to form the mass matrix and
    !! the external force vector.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates at which to evaluate the matrix.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 6-by-6N shape function matrix, where N is the number of nodes.

    ! Local Variables
    integer(int32) :: i, k, n
    real(real64) :: ni

    ! Process
    n = this%get_node_count()
    allocate(rst(6, 6 * n), source = 0.0d0)
    do i = 1, n
        ni = this%evaluate_shape_function(i, s)
        do k = 1, 6
            rst(k, 6 * (i - 1) + k) = ni
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function shl_area(this) result(rst)
    !! Computes the area of the element.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64) :: rst
        !! The element area.

    ! Local Variables
    integer(int32) :: i
    real(real64), allocatable, dimension(:) :: wts
    real(real64), allocatable, dimension(:,:) :: pts

    ! Process
    call this%integration_rule(pts, wts)
    rst = 0.0d0
    do i = 1, size(wts)
        rst = rst + wts(i) * det2(this%jacobian(pts(:,i)))
    end do
end function

! ------------------------------------------------------------------------------
pure function shl_stiffness_matrix(this, rule) result(rst)
    !! Computes the element stiffness matrix in the global coordinate system.
    !! The local stiffness matrix is
    !! $$ K_{local} = \int_A B^T D B \, dA + \int_A \alpha G t \,
    !! \boldsymbol{b}_d \boldsymbol{b}_d^T dA, $$
    !! where \(\boldsymbol{b}_d\) is the drilling-rotation operator.  The
    !! global matrix is \(K = T^T K_{local} T\).
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused and is present for
        !! interface compatibility; the element uses the quadrature defined by
        !! its integration_rule routine.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 6N-by-6N stiffness matrix, where N is the number of nodes.

    ! Local Variables
    integer(int32) :: i, ndof
    real(real64) :: dj, g, kd
    real(real64), allocatable, dimension(:) :: wts, bd
    real(real64), allocatable, dimension(:,:) :: pts, b, d, t, kl

    ! Initialization
    call this%integration_rule(pts, wts)
    d = this%constitutive_matrix()
    ndof = 6 * this%get_node_count()
    g = this%material%modulus / (2.0d0 * (1.0d0 + &
        this%material%poissons_ratio))
    kd = this%drilling_factor * g * this%thickness
    allocate(kl(ndof, ndof), source = 0.0d0)

    ! Integrate the membrane, bending, shear, and drilling contributions
    do i = 1, size(wts)
        dj = wts(i) * det2(this%jacobian(pts(:,i)))
        b = this%strain_displacement_matrix(pts(:,i))
        kl = kl + dj * matmul(transpose(b), matmul(d, b))
        bd = drilling_row(this, pts(:,i))
        kl = kl + (dj * kd) * spread(bd, 2, ndof) * spread(bd, 1, ndof)
    end do

    ! Transform to global coordinates, removing round-off asymmetry
    t = this%rotation_matrix()
    rst = matmul(transpose(t), matmul(kl, t))
    rst = 0.5d0 * (rst + transpose(rst))
end function

! ------------------------------------------------------------------------------
pure function shl_mass_matrix(this, rule) result(rst)
    !! Computes the consistent element mass matrix in the global coordinate
    !! system.
    !! $$ M = T^T \left( \int_A \rho N^T \Lambda N \, dA \right) T, \quad
    !! \Lambda = \mathrm{diag}\left(t, t, t, \frac{t^3}{12},
    !! \frac{t^3}{12}, \frac{t^3}{12}\right) $$
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused and is present for
        !! interface compatibility; the element uses the quadrature defined by
        !! its integration_rule routine.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 6N-by-6N mass matrix, where N is the number of nodes.

    ! Local Variables
    integer(int32) :: i, k, ndof
    real(real64) :: dj, lambda(6)
    real(real64), allocatable, dimension(:) :: wts
    real(real64), allocatable, dimension(:,:) :: pts, n, ln, t, ml

    ! Initialization
    call this%integration_rule(pts, wts)
    ndof = 6 * this%get_node_count()
    lambda(1:3) = this%thickness
    lambda(4:6) = this%thickness**3 / 12.0d0
    lambda = this%material%density * lambda
    allocate(ml(ndof, ndof), source = 0.0d0)

    ! Integrate
    do i = 1, size(wts)
        dj = wts(i) * det2(this%jacobian(pts(:,i)))
        n = this%shape_function_matrix(pts(:,i))
        ln = n
        do k = 1, 6
            ln(k,:) = lambda(k) * n(k,:)
        end do
        ml = ml + dj * matmul(transpose(n), ln)
    end do

    ! Transform to global coordinates, removing round-off asymmetry
    t = this%rotation_matrix()
    rst = matmul(transpose(t), matmul(ml, t))
    rst = 0.5d0 * (rst + transpose(rst))
end function

! ------------------------------------------------------------------------------
pure function shl_ext_force_vector(this, q, rule) result(rst)
    !! Computes the consistent nodal force vector, in the global coordinate
    !! system, resulting from a uniform surface traction.
    !! $$ \boldsymbol{f}_i = \left( \int_A N_i \, dA \right) \boldsymbol{q} $$
    !! Only the translational degrees of freedom receive load.
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), intent(in), dimension(:) :: q
        !! The 3-element surface traction vector (force per unit area)
        !! expressed in the global coordinate system.  For a pressure \(p\)
        !! acting along the local \(z\)-axis, supply
        !! \(\boldsymbol{q} = p \boldsymbol{e}_z\), where
        !! \(\boldsymbol{e}_z\) is the third row of the local_frame matrix.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  This argument is unused and is present for
        !! interface compatibility; the element uses the quadrature defined by
        !! its integration_rule routine.
    real(real64), allocatable, dimension(:) :: rst
        !! The 6N-element force vector, where N is the number of nodes.

    ! Local Variables
    integer(int32) :: i, j, c, n
    real(real64) :: dj
    real(real64), allocatable, dimension(:) :: wts
    real(real64), allocatable, dimension(:,:) :: pts

    ! Input Checking
    if (size(q) /= 3) error stop DYN_ARRAY_SIZE_ERROR

    ! Process
    call this%integration_rule(pts, wts)
    n = this%get_node_count()
    allocate(rst(6 * n), source = 0.0d0)
    do i = 1, size(wts)
        dj = wts(i) * det2(this%jacobian(pts(:,i)))
        do j = 1, n
            c = 6 * (j - 1)
            rst(c+1:c+3) = rst(c+1:c+3) + &
                dj * this%evaluate_shape_function(j, pts(:,i)) * q
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function shl_strain(this, displacement, s) result(rst)
    !! Computes the generalized strain vector, in the local coordinate
    !! system, at the specified natural coordinate.  The components are
    !! ordered as
    !! \([\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}, \kappa_{xx},
    !! \kappa_{yy}, \kappa_{xy}, \gamma_{xz}, \gamma_{yz}]\).
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), intent(in), dimension(:) :: displacement
        !! The 6N-element displacement vector in the global coordinate
        !! system, where N is the number of nodes.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates at which to evaluate the strain.
    real(real64), allocatable, dimension(:) :: rst
        !! The 8-element generalized strain vector.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: t

    ! Process
    t = this%rotation_matrix()
    if (size(displacement) /= size(t, 2)) error stop DYN_ARRAY_SIZE_ERROR
    rst = matmul(this%strain_displacement_matrix(s), &
        matmul(t, displacement))
end function

! ------------------------------------------------------------------------------
pure function shl_stress(this, displacement, s) result(rst)
    !! Computes the stress-resultant vector, in the local coordinate system,
    !! at the specified natural coordinate.  The components are ordered as
    !! \([N_{xx}, N_{yy}, N_{xy}, M_{xx}, M_{yy}, M_{xy}, Q_x, Q_y]\).  The
    !! surface stresses may be recovered as
    !! \(\sigma = N / t \pm 6 M / t^2\).
    class(shell_element), intent(in) :: this
        !! The shell_element object.
    real(real64), intent(in), dimension(:) :: displacement
        !! The 6N-element displacement vector in the global coordinate
        !! system, where N is the number of nodes.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates at which to evaluate the stress resultants.
    real(real64), allocatable, dimension(:) :: rst
        !! The 8-element stress-resultant vector.
    rst = matmul(this%constitutive_matrix(), this%strain(displacement, s))
end function

! ******************************************************************************
! TRIANGULAR_SHELL_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function tri_init(mat, thickness, nd1, nd2, nd3) result(rst)
    !! Constructs a new [[triangular_shell_element]].
    class(material), intent(in) :: mat
        !! The material.
    real(real64), intent(in) :: thickness
        !! The shell thickness.  This value must be positive.
    class(node), intent(in) :: nd1
        !! The first node.
    class(node), intent(in) :: nd2
        !! The second node.
    class(node), intent(in) :: nd3
        !! The third node.
    type(triangular_shell_element) :: rst
        !! The new [[triangular_shell_element]].

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: frame

    ! Input Checking
    if (thickness <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR

    ! Process
    rst%material = mat
    rst%thickness = thickness
    rst%nodes(1) = nd1
    rst%nodes(2) = nd2
    rst%nodes(3) = nd3

    ! Verify the geometry is not degenerate
    frame = rst%local_frame()
end function

! ------------------------------------------------------------------------------
pure function tri_get_node_count(this) result(rst)
    !! Gets the number of nodes in the element.
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    integer(int32) :: rst
        !! The number of nodes (always 3).
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function tri_get_node(this, i) result(rst)
    !! Gets the requested node from the element.
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node to retrieve.
    type(node) :: rst
        !! The requested node.
    if (i < 1 .or. i > 3) error stop DYN_INDEX_OUT_OF_RANGE
    rst = this%nodes(i)
end function

! ------------------------------------------------------------------------------
pure function tri_get_node_natural_coordinates(this, i) result(rst)
    !! Returns the natural coordinates of the requested node.
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node.
    real(real64), allocatable, dimension(:) :: rst
        !! The 2-element natural coordinate vector of the node.
    select case (i)
    case (1)
        rst = [0.0d0, 0.0d0]
    case (2)
        rst = [1.0d0, 0.0d0]
    case (3)
        rst = [0.0d0, 1.0d0]
    case default
        error stop DYN_INDEX_OUT_OF_RANGE
    end select
end function

! ------------------------------------------------------------------------------
pure function tri_shape_function(this, i, s) result(rst)
    !! Evaluates the i-th linear shape function,
    !! \(N_1 = 1 - \xi - \eta\), \(N_2 = \xi\), \(N_3 = \eta\).
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    integer(int32), intent(in) :: i
        !! The index of the shape function to evaluate.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates \((\xi, \eta)\).
    real(real64) :: rst
        !! The value of the i-th shape function.
    select case (i)
    case (1)
        rst = 1.0d0 - s(1) - s(2)
    case (2)
        rst = s(1)
    case (3)
        rst = s(2)
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure function tri_shape_function_natural_gradient(this, s) result(rst)
    !! Computes the derivatives of the linear shape functions with respect
    !! to the natural coordinates.
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates \((\xi, \eta)\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 2-by-3 matrix of derivatives.
    rst = reshape([-1.0d0, -1.0d0, 1.0d0, 0.0d0, 0.0d0, 1.0d0], [2, 3])
end function

! ------------------------------------------------------------------------------
pure subroutine tri_integration_rule(this, pts, wts)
    !! Returns the three-point interior quadrature rule for triangles, which
    !! integrates quadratic polynomials exactly.
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    real(real64), allocatable, intent(out), dimension(:,:) :: pts
        !! The 2-by-3 matrix of integration point natural coordinates.
    real(real64), allocatable, intent(out), dimension(:) :: wts
        !! The 3 integration weights.

    real(real64), parameter :: a = 1.0d0 / 6.0d0
    real(real64), parameter :: b = 2.0d0 / 3.0d0

    pts = reshape([a, a, b, a, a, b], [2, 3])
    wts = [a, a, a]
end subroutine

! ------------------------------------------------------------------------------
pure function dkt_curvature_matrix(xy, s) result(rst)
    ! Computes the 3-by-9 DKT curvature-displacement matrix, relating the
    ! nodal [w, theta_x, theta_y] values to [kxx, kyy, kxy], following Batoz,
    ! Bathe, and Ho (1980).  The DKT rotations are beta_x = theta_y and
    ! beta_y = -theta_x.
    real(real64), intent(in), dimension(3,2) :: xy
        ! The local nodal coordinates.
    real(real64), intent(in), dimension(:) :: s
        ! The natural coordinates (xi, eta).
    real(real64) :: rst(3,9)

    ! Edge k runs from node ii(k) to node jj(k); edges 1, 2, and 3 correspond
    ! to the midside nodes 4 (2-3), 5 (3-1), and 6 (1-2) of Batoz et al.
    integer(int32), parameter :: ii(3) = [2, 3, 1]
    integer(int32), parameter :: jj(3) = [3, 1, 2]

    ! Local Variables
    integer(int32) :: k, m
    real(real64) :: xi, eta, xij, yij, l2, x31, x12, y31, y12, twoa
    real(real64) :: a(3), b(3), c(3), d(3), e(3), dn(2,6), hx(2,9), hy(2,9)

    ! Edge geometry coefficients
    do k = 1, 3
        xij = xy(ii(k),1) - xy(jj(k),1)
        yij = xy(ii(k),2) - xy(jj(k),2)
        l2 = xij**2 + yij**2
        a(k) = -xij / l2
        b(k) = 0.75d0 * xij * yij / l2
        c(k) = (0.25d0 * xij**2 - 0.5d0 * yij**2) / l2
        d(k) = -yij / l2
        e(k) = (0.25d0 * yij**2 - 0.5d0 * xij**2) / l2
    end do

    ! Derivatives of the 6-node quadratic shape functions with respect to xi
    ! (row 1) and eta (row 2)
    xi = s(1)
    eta = s(2)
    dn(1,:) = [4.0d0 * (xi + eta) - 3.0d0, 4.0d0 * xi - 1.0d0, 0.0d0, &
        4.0d0 * eta, -4.0d0 * eta, 4.0d0 * (1.0d0 - 2.0d0 * xi - eta)]
    dn(2,:) = [4.0d0 * (xi + eta) - 3.0d0, 0.0d0, 4.0d0 * eta - 1.0d0, &
        4.0d0 * xi, 4.0d0 * (1.0d0 - xi - 2.0d0 * eta), -4.0d0 * xi]

    ! Derivatives of the beta_x (Hx) and beta_y (Hy) interpolation functions
    do m = 1, 2
        hx(m,1) = 1.5d0 * (a(3) * dn(m,6) - a(2) * dn(m,5))
        hx(m,2) = b(2) * dn(m,5) + b(3) * dn(m,6)
        hx(m,3) = dn(m,1) - c(2) * dn(m,5) - c(3) * dn(m,6)
        hx(m,4) = 1.5d0 * (a(1) * dn(m,4) - a(3) * dn(m,6))
        hx(m,5) = b(3) * dn(m,6) + b(1) * dn(m,4)
        hx(m,6) = dn(m,2) - c(3) * dn(m,6) - c(1) * dn(m,4)
        hx(m,7) = 1.5d0 * (a(2) * dn(m,5) - a(1) * dn(m,4))
        hx(m,8) = b(1) * dn(m,4) + b(2) * dn(m,5)
        hx(m,9) = dn(m,3) - c(1) * dn(m,4) - c(2) * dn(m,5)

        hy(m,1) = 1.5d0 * (d(3) * dn(m,6) - d(2) * dn(m,5))
        hy(m,2) = -dn(m,1) + e(2) * dn(m,5) + e(3) * dn(m,6)
        hy(m,3) = -hx(m,2)
        hy(m,4) = 1.5d0 * (d(1) * dn(m,4) - d(3) * dn(m,6))
        hy(m,5) = -dn(m,2) + e(3) * dn(m,6) + e(1) * dn(m,4)
        hy(m,6) = -hx(m,5)
        hy(m,7) = 1.5d0 * (d(2) * dn(m,5) - d(1) * dn(m,4))
        hy(m,8) = -dn(m,3) + e(1) * dn(m,4) + e(2) * dn(m,5)
        hy(m,9) = -hx(m,8)
    end do

    ! Map natural derivatives to local x-y derivatives
    x31 = xy(3,1) - xy(1,1)
    x12 = xy(1,1) - xy(2,1)
    y31 = xy(3,2) - xy(1,2)
    y12 = xy(1,2) - xy(2,2)
    twoa = x31 * y12 - x12 * y31
    rst(1,:) = (y31 * hx(1,:) + y12 * hx(2,:)) / twoa
    rst(2,:) = (-x31 * hy(1,:) - x12 * hy(2,:)) / twoa
    rst(3,:) = (-x31 * hx(1,:) - x12 * hx(2,:) + y31 * hy(1,:) + &
        y12 * hy(2,:)) / twoa
end function

! ------------------------------------------------------------------------------
pure function tri_strain_disp_matrix(this, s) result(rst)
    !! Computes the 8-by-18 generalized strain-displacement matrix in the
    !! local coordinate system.  Rows 1-3 contain the constant-strain
    !! membrane terms, rows 4-6 the DKT curvature terms, and rows 7-8 (the
    !! transverse shear strains) are zero as the DKT formulation neglects
    !! transverse shear deformation.
    class(triangular_shell_element), intent(in) :: this
        !! The triangular_shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates \((\xi, \eta)\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 8-by-18 strain-displacement matrix.

    ! Local Variables
    integer(int32) :: i, m
    real(real64) :: bb(3,9)
    real(real64), allocatable, dimension(:,:) :: xy

    ! Membrane terms
    allocate(rst(8, 18), source = 0.0d0)
    rst(1:3,:) = membrane_strain_rows(this%shape_function_gradient(s))

    ! Bending terms: DKT local DOF (w, theta_x, theta_y) of node i map onto
    ! element DOF 3, 4, and 5 of the node.
    xy = this%local_coordinates()
    bb = dkt_curvature_matrix(xy, s)
    do i = 1, 3
        do m = 1, 3
            rst(4:6, 6 * (i - 1) + 2 + m) = bb(:, 3 * (i - 1) + m)
        end do
    end do
end function

! ******************************************************************************
! RECTANGULAR_SHELL_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function quad_init(mat, thickness, nd1, nd2, nd3, nd4) result(rst)
    !! Constructs a new [[rectangular_shell_element]].
    class(material), intent(in) :: mat
        !! The material.
    real(real64), intent(in) :: thickness
        !! The shell thickness.  This value must be positive.
    class(node), intent(in) :: nd1
        !! The first node.
    class(node), intent(in) :: nd2
        !! The second node.
    class(node), intent(in) :: nd3
        !! The third node.
    class(node), intent(in) :: nd4
        !! The fourth node.
    type(rectangular_shell_element) :: rst
        !! The new [[rectangular_shell_element]].

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: frame

    ! Input Checking
    if (thickness <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR

    ! Process
    rst%material = mat
    rst%thickness = thickness
    rst%nodes(1) = nd1
    rst%nodes(2) = nd2
    rst%nodes(3) = nd3
    rst%nodes(4) = nd4

    ! Verify the geometry is not degenerate
    frame = rst%local_frame()
end function

! ------------------------------------------------------------------------------
pure function quad_get_node_count(this) result(rst)
    !! Gets the number of nodes in the element.
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    integer(int32) :: rst
        !! The number of nodes (always 4).
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function quad_get_node(this, i) result(rst)
    !! Gets the requested node from the element.
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node to retrieve.
    type(node) :: rst
        !! The requested node.
    if (i < 1 .or. i > 4) error stop DYN_INDEX_OUT_OF_RANGE
    rst = this%nodes(i)
end function

! ------------------------------------------------------------------------------
pure function quad_get_node_natural_coordinates(this, i) result(rst)
    !! Returns the natural coordinates of the requested node.
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    integer(int32), intent(in) :: i
        !! The local index of the node.
    real(real64), allocatable, dimension(:) :: rst
        !! The 2-element natural coordinate vector of the node.

    real(real64), parameter :: rn(4) = [-1.0d0, 1.0d0, 1.0d0, -1.0d0]
    real(real64), parameter :: sn(4) = [-1.0d0, -1.0d0, 1.0d0, 1.0d0]

    if (i < 1 .or. i > 4) error stop DYN_INDEX_OUT_OF_RANGE
    rst = [rn(i), sn(i)]
end function

! ------------------------------------------------------------------------------
pure function quad_shape_function(this, i, s) result(rst)
    !! Evaluates the i-th bilinear shape function,
    !! \(N_i = \frac{1}{4}(1 + r r_i)(1 + s s_i)\).
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    integer(int32), intent(in) :: i
        !! The index of the shape function to evaluate.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates \((r, s)\).
    real(real64) :: rst
        !! The value of the i-th shape function.

    real(real64), parameter :: rn(4) = [-1.0d0, 1.0d0, 1.0d0, -1.0d0]
    real(real64), parameter :: sn(4) = [-1.0d0, -1.0d0, 1.0d0, 1.0d0]

    if (i < 1 .or. i > 4) then
        rst = 0.0d0
    else
        rst = 0.25d0 * (1.0d0 + s(1) * rn(i)) * (1.0d0 + s(2) * sn(i))
    end if
end function

! ------------------------------------------------------------------------------
pure function quad_shape_function_natural_gradient(this, s) result(rst)
    !! Computes the derivatives of the bilinear shape functions with respect
    !! to the natural coordinates.
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates \((r, s)\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 2-by-4 matrix of derivatives.

    real(real64), parameter :: rn(4) = [-1.0d0, 1.0d0, 1.0d0, -1.0d0]
    real(real64), parameter :: sn(4) = [-1.0d0, -1.0d0, 1.0d0, 1.0d0]

    allocate(rst(2, 4))
    rst(1,:) = 0.25d0 * rn * (1.0d0 + s(2) * sn)
    rst(2,:) = 0.25d0 * sn * (1.0d0 + s(1) * rn)
end function

! ------------------------------------------------------------------------------
pure subroutine quad_integration_rule(this, pts, wts)
    !! Returns the 2-by-2 Gauss quadrature rule.
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    real(real64), allocatable, intent(out), dimension(:,:) :: pts
        !! The 2-by-4 matrix of integration point natural coordinates.
    real(real64), allocatable, intent(out), dimension(:) :: wts
        !! The 4 integration weights.

    real(real64) :: g

    g = 1.0d0 / sqrt(3.0d0)
    pts = reshape([-g, -g, g, -g, g, g, -g, g], [2, 4])
    wts = [1.0d0, 1.0d0, 1.0d0, 1.0d0]
end subroutine

! ------------------------------------------------------------------------------
pure function mitc4_covariant_shear(elem, r, s, dir) result(rst)
    ! Computes the 24-element operator giving the covariant transverse shear
    ! strain at (r, s):
    !   dir = 1: gamma_r = dw/dr + theta_y * dx/dr - theta_x * dy/dr
    !   dir = 2: gamma_s = dw/ds + theta_y * dx/ds - theta_x * dy/ds
    class(rectangular_shell_element), intent(in) :: elem
    real(real64), intent(in) :: r, s
    integer(int32), intent(in) :: dir
    real(real64) :: rst(24)

    integer(int32) :: i, c
    real(real64) :: ni, jac(2,2)
    real(real64), allocatable, dimension(:,:) :: dnds

    dnds = elem%shape_function_natural_gradient([r, s])
    jac = elem%jacobian([r, s])
    rst = 0.0d0
    do i = 1, 4
        c = 6 * (i - 1)
        ni = elem%evaluate_shape_function(i, [r, s])
        rst(c+3) = dnds(dir,i)
        rst(c+4) = -ni * jac(dir,2)
        rst(c+5) = ni * jac(dir,1)
    end do
end function

! ------------------------------------------------------------------------------
pure function quad_strain_disp_matrix(this, s) result(rst)
    !! Computes the 8-by-24 generalized strain-displacement matrix in the
    !! local coordinate system.  Rows 1-3 contain the bilinear membrane
    !! terms, rows 4-6 the Mindlin-Reissner curvature terms, and rows 7-8 the
    !! MITC4 assumed transverse shear strains.
    !!
    !! The covariant transverse shear strains are sampled at the edge
    !! midpoints \(A = (0, 1)\), \(B = (-1, 0)\), \(C = (0, -1)\), and
    !! \(D = (1, 0)\) and interpolated as
    !! $$ \tilde{\gamma}_r = \frac{1}{2}(1 + s)\gamma_r^A +
    !! \frac{1}{2}(1 - s)\gamma_r^C, \quad
    !! \tilde{\gamma}_s = \frac{1}{2}(1 + r)\gamma_s^D +
    !! \frac{1}{2}(1 - r)\gamma_s^B, $$
    !! and then transformed to local Cartesian components with
    !! \([\gamma_{xz}, \gamma_{yz}]^T = J^{-1}[\tilde{\gamma}_r,
    !! \tilde{\gamma}_s]^T\).
    class(rectangular_shell_element), intent(in) :: this
        !! The rectangular_shell_element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinates \((r, s)\).
    real(real64), allocatable, dimension(:,:) :: rst
        !! The 8-by-24 strain-displacement matrix.

    ! Local Variables
    integer(int32) :: i, c
    real(real64) :: jac(2,2), jinv(2,2), gr(24), gs(24)
    real(real64), allocatable, dimension(:,:) :: dndx

    ! Membrane terms
    allocate(rst(8, 24), source = 0.0d0)
    dndx = this%shape_function_gradient(s)
    rst(1:3,:) = membrane_strain_rows(dndx)

    ! Bending terms
    do i = 1, 4
        c = 6 * (i - 1)
        rst(4,c+5) = dndx(1,i)
        rst(5,c+4) = -dndx(2,i)
        rst(6,c+4) = -dndx(1,i)
        rst(6,c+5) = dndx(2,i)
    end do

    ! Assumed covariant transverse shear strains
    gr = 0.5d0 * (1.0d0 + s(2)) * &
        mitc4_covariant_shear(this, 0.0d0, 1.0d0, 1) + &
        0.5d0 * (1.0d0 - s(2)) * &
        mitc4_covariant_shear(this, 0.0d0, -1.0d0, 1)
    gs = 0.5d0 * (1.0d0 + s(1)) * &
        mitc4_covariant_shear(this, 1.0d0, 0.0d0, 2) + &
        0.5d0 * (1.0d0 - s(1)) * &
        mitc4_covariant_shear(this, -1.0d0, 0.0d0, 2)

    ! Transform to local Cartesian components
    jac = this%jacobian(s)
    jinv = reshape([jac(2,2), -jac(2,1), -jac(1,2), jac(1,1)], [2, 2]) / &
        det2(jac)
    rst(7,:) = jinv(1,1) * gr + jinv(1,2) * gs
    rst(8,:) = jinv(2,1) * gr + jinv(2,2) * gs
end function

! ------------------------------------------------------------------------------
end module
