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
module dynamics_truss_elements
    !! Two-node, axial-only planar and spatial truss elements.
    !! Each node carries only translational degrees of freedom. The global
    !! stiffness of a member of length L and unit direction n is
    !!
    !! $$K = \frac{E A}{L}\begin{bmatrix}-n\\n\end{bmatrix}
    !! \begin{bmatrix}-n^T & n^T\end{bmatrix}.$$
    !!
    !! The inherited mass matrix uses consistent linear interpolation,
    !! giving rho*A*L/6 times [[2I,I],[I,2I]]. Neither element resists
    !! bending or has nodal rotational degrees of freedom.
    use iso_fortran_env, only : int32, real64
    use dynamics_structural, only : material, node, line_element
    use dynamics_helper, only : cross_product
    implicit none
    private
    public :: truss_element_2d, truss_element_3d

! ------------------------------------------------------------------------------
    type, extends(line_element) :: truss_element_2d
        !! Defines a pin-jointed 2D bar with x and y translations per node.
        !! Use the parent mass_matrix and the global axial stiffness override.
        type(node) :: node_1
            !! The node at natural coordinate s = -1.
        type(node) :: node_2
            !! The node at natural coordinate s = 1.
    contains
        procedure, public :: get_dimensionality => t2d_dimensionality
        procedure, public :: get_node_count => t2d_node_count
        procedure, public :: get_dof_per_node => t2d_dof_per_node
        procedure, public :: get_node => t2d_get_node
        procedure, public :: get_terminal_nodes => t2d_terminal_nodes
        procedure, public :: evaluate_shape_function => t2d_shape_function
        procedure, public :: shape_function_matrix => t2d_shape_matrix
        procedure, public :: strain_displacement_matrix => t2d_strain_matrix
        procedure, public :: constitutive_matrix => t2d_constitutive_matrix
        procedure, public :: jacobian => t2d_jacobian
        procedure, public :: rotation_matrix => t2d_rotation_matrix
        procedure, public :: stiffness_matrix => t2d_stiffness_matrix
    end type

    interface truss_element_2d
        module procedure :: t2d_init
    end interface

! ------------------------------------------------------------------------------
    type, extends(line_element) :: truss_element_3d
        !! Defines a pin-jointed 3D bar with three translations per node.
        !! The local transverse axes are arbitrary; only the longitudinal
        !! direction influences axial strain and stiffness.
        type(node) :: node_1
            !! The node at natural coordinate s = -1.
        type(node) :: node_2
            !! The node at natural coordinate s = 1.
    contains
        procedure, public :: get_dimensionality => t3d_dimensionality
        procedure, public :: get_node_count => t3d_node_count
        procedure, public :: get_dof_per_node => t3d_dof_per_node
        procedure, public :: get_node => t3d_get_node
        procedure, public :: get_terminal_nodes => t3d_terminal_nodes
        procedure, public :: evaluate_shape_function => t3d_shape_function
        procedure, public :: shape_function_matrix => t3d_shape_matrix
        procedure, public :: strain_displacement_matrix => t3d_strain_matrix
        procedure, public :: constitutive_matrix => t3d_constitutive_matrix
        procedure, public :: jacobian => t3d_jacobian
        procedure, public :: rotation_matrix => t3d_rotation_matrix
        procedure, public :: stiffness_matrix => t3d_stiffness_matrix
    end type

    interface truss_element_3d
        module procedure :: t3d_init
    end interface

contains
! ******************************************************************************
! 2D TRUSS ELEMENT
! ------------------------------------------------------------------------------
pure function t2d_init(mat, area, nd1, nd2) result(rst)
    !! Initializes a [[truss_element_2d]].
    class(material), intent(in) :: mat
    !! The elastic material and density.
    real(real64), intent(in) :: area
    !! The cross-sectional area.
    class(node), intent(in) :: nd1, nd2
    !! The first and second truss nodes.
    type(truss_element_2d) :: rst
    !! The initialized planar truss element.
    rst%material = mat
    rst%area = area
    rst%node_1 = nd1
    rst%node_2 = nd2
end function

! ------------------------------------------------------------------------------
pure function t2d_dimensionality(this) result(rst)
    class(truss_element_2d), intent(in) :: this
    integer(int32) :: rst
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function t2d_node_count(this) result(rst)
    class(truss_element_2d), intent(in) :: this
    integer(int32) :: rst
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function t2d_dof_per_node(this) result(rst)
    class(truss_element_2d), intent(in) :: this
    integer(int32) :: rst
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function t2d_get_node(this, i) result(rst)
    class(truss_element_2d), intent(in) :: this
    integer(int32), intent(in) :: i
    type(node) :: rst
    if (i == 1) then
        rst = this%node_1
    else
        rst = this%node_2
    end if
end function

! ------------------------------------------------------------------------------
pure subroutine t2d_terminal_nodes(this, i1, i2)
    class(truss_element_2d), intent(in) :: this
    integer(int32), intent(out) :: i1, i2
    i1 = 1
    i2 = 2
end subroutine

! ------------------------------------------------------------------------------
pure function t2d_shape_function(this, i, s) result(rst)
    !! Linear interpolation in the natural coordinate -1 <= s <= 1.
    class(truss_element_2d), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64), intent(in), dimension(:) :: s
    real(real64) :: rst
    select case (i)
    case (1)
        rst = 0.5d0 * (1.0d0 - s(1))
    case (2)
        rst = 0.5d0 * (1.0d0 + s(1))
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure function t2d_shape_matrix(this, s) result(rst)
    !! Computes the two-component linear displacement interpolation matrix.
    class(truss_element_2d), intent(in) :: this
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: rst
    real(real64) :: n1, n2
    n1 = this%evaluate_shape_function(1, s)
    n2 = this%evaluate_shape_function(2, s)
    allocate(rst(2,4), source = 0.0d0)
    rst(1,1) = n1
    rst(2,2) = n1
    rst(1,3) = n2
    rst(2,4) = n2
end function

! ------------------------------------------------------------------------------
pure function t2d_strain_matrix(this, s) result(rst)
    !! Local axial strain is (u_2 - u_1) / L.
    class(truss_element_2d), intent(in) :: this
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: rst
    allocate(rst(1,4), source = 0.0d0)
    rst(1,1) = -1.0d0 / this%length()
    rst(1,3) = -rst(1,1)
end function

! ------------------------------------------------------------------------------
pure function t2d_constitutive_matrix(this) result(rst)
    !! Maps axial strain to axial force (E A times strain).
    class(truss_element_2d), intent(in) :: this
    real(real64), allocatable, dimension(:,:) :: rst
    allocate(rst(1,1), source = this%material%modulus * this%area)
end function

! ------------------------------------------------------------------------------
pure function t2d_jacobian(this, s) result(rst)
    !! The natural-to-physical coordinate Jacobian is L/2.
    class(truss_element_2d), intent(in) :: this
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: rst
    allocate(rst(1,1), source = 0.5d0 * this%length())
end function

! ------------------------------------------------------------------------------
pure function t2d_rotation_matrix(this) result(rst)
    !! Maps local translations to global translations at each node.
    class(truss_element_2d), intent(in) :: this
    real(real64), allocatable, dimension(:,:) :: rst
    real(real64) :: cosine, sine, length
    length = this%length()
    cosine = (this%node_2%x - this%node_1%x) / length
    sine = (this%node_2%y - this%node_1%y) / length
    allocate(rst(4,4), source = 0.0d0)
    rst(1:2,1:2) = reshape([cosine, sine, -sine, cosine], [2,2])
    rst(3:4,3:4) = rst(1:2,1:2)
end function

! ------------------------------------------------------------------------------
pure function t2d_stiffness_matrix(this, rule) result(rst)
    !! Global axial stiffness K = (E A/L) b b^T, where b = [-n, n].
    class(truss_element_2d), intent(in) :: this
    !! The planar truss element.
    integer(int32), intent(in), optional :: rule
    !! The integration rule is accepted for compatibility; stiffness is exact.
    real(real64), allocatable, dimension(:,:) :: rst
    !! The 4-by-4 global translational stiffness matrix.
    real(real64) :: length, direction(2), axial(4)
    integer(int32) :: row, col
    length = this%length()
    direction = [this%node_2%x - this%node_1%x, &
        this%node_2%y - this%node_1%y] / length
    axial = [-direction, direction]
    allocate(rst(4,4))
    do col = 1, 4
        do row = 1, 4
            rst(row,col) = this%material%modulus * this%area / length * axial(row) * axial(col)
        end do
    end do
end function

! ******************************************************************************
! 3D TRUSS ELEMENT
! ------------------------------------------------------------------------------
pure function t3d_init(mat, area, nd1, nd2) result(rst)
    !! Initializes a [[truss_element_3d]].
    class(material), intent(in) :: mat
    !! The elastic material and density.
    real(real64), intent(in) :: area
    !! The cross-sectional area.
    class(node), intent(in) :: nd1, nd2
    !! The first and second truss nodes.
    type(truss_element_3d) :: rst
    !! The initialized spatial truss element.
    rst%material = mat
    rst%area = area
    rst%node_1 = nd1
    rst%node_2 = nd2
end function

! ------------------------------------------------------------------------------
pure function t3d_dimensionality(this) result(rst)
    class(truss_element_3d), intent(in) :: this
    integer(int32) :: rst
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function t3d_node_count(this) result(rst)
    class(truss_element_3d), intent(in) :: this
    integer(int32) :: rst
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function t3d_dof_per_node(this) result(rst)
    class(truss_element_3d), intent(in) :: this
    integer(int32) :: rst
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function t3d_get_node(this, i) result(rst)
    class(truss_element_3d), intent(in) :: this
    integer(int32), intent(in) :: i
    type(node) :: rst
    if (i == 1) then
        rst = this%node_1
    else
        rst = this%node_2
    end if
end function

! ------------------------------------------------------------------------------
pure subroutine t3d_terminal_nodes(this, i1, i2)
    class(truss_element_3d), intent(in) :: this
    integer(int32), intent(out) :: i1, i2
    i1 = 1
    i2 = 2
end subroutine

! ------------------------------------------------------------------------------
pure function t3d_shape_function(this, i, s) result(rst)
    class(truss_element_3d), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64), intent(in), dimension(:) :: s
    real(real64) :: rst
    select case (i)
    case (1)
        rst = 0.5d0 * (1.0d0 - s(1))
    case (2)
        rst = 0.5d0 * (1.0d0 + s(1))
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure function t3d_shape_matrix(this, s) result(rst)
    !! Computes the three-component linear displacement interpolation matrix.
    class(truss_element_3d), intent(in) :: this
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: rst
    real(real64) :: n1, n2
    integer(int32) :: axis
    n1 = this%evaluate_shape_function(1, s)
    n2 = this%evaluate_shape_function(2, s)
    allocate(rst(3,6), source = 0.0d0)
    do axis = 1, 3
        rst(axis,axis) = n1
        rst(axis,axis+3) = n2
    end do
end function

! ------------------------------------------------------------------------------
pure function t3d_strain_matrix(this, s) result(rst)
    !! Local axial strain is (u_2 - u_1) / L.
    class(truss_element_3d), intent(in) :: this
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: rst
    allocate(rst(1,6), source = 0.0d0)
    rst(1,1) = -1.0d0 / this%length()
    rst(1,4) = -rst(1,1)
end function

! ------------------------------------------------------------------------------
pure function t3d_constitutive_matrix(this) result(rst)
    !! Maps axial strain to axial force (E A times strain).
    class(truss_element_3d), intent(in) :: this
    real(real64), allocatable, dimension(:,:) :: rst
    allocate(rst(1,1), source = this%material%modulus * this%area)
end function

! ------------------------------------------------------------------------------
pure function t3d_jacobian(this, s) result(rst)
    !! The natural-to-physical coordinate Jacobian is L/2.
    class(truss_element_3d), intent(in) :: this
    real(real64), intent(in), dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: rst
    allocate(rst(1,1), source = 0.5d0 * this%length())
end function

! ------------------------------------------------------------------------------
pure function t3d_rotation_matrix(this) result(rst)
    !! Returns a local-to-global orthonormal basis with the first axis along
    !! the bar. Transverse axes do not affect its axial stiffness.
    class(truss_element_3d), intent(in) :: this
    real(real64), allocatable, dimension(:,:) :: rst
    real(real64) :: longitudinal(3), transverse(3), normal(3), reference(3)
    integer(int32) :: end_index
    longitudinal = [this%node_2%x - this%node_1%x, &
        this%node_2%y - this%node_1%y, this%node_2%z - this%node_1%z] / this%length()
    if (abs(longitudinal(3)) < 0.9d0) then
        reference = [0.0d0, 0.0d0, 1.0d0]
    else
        reference = [0.0d0, 1.0d0, 0.0d0]
    end if
    transverse = cross_product(reference, longitudinal)
    transverse = transverse / norm2(transverse)
    normal = cross_product(longitudinal, transverse)
    allocate(rst(6,6), source = 0.0d0)
    do end_index = 0, 1
        rst(1+3*end_index:3+3*end_index,1+3*end_index) = longitudinal
        rst(1+3*end_index:3+3*end_index,2+3*end_index) = transverse
        rst(1+3*end_index:3+3*end_index,3+3*end_index) = normal
    end do
end function

! ------------------------------------------------------------------------------
pure function t3d_stiffness_matrix(this, rule) result(rst)
    !! Global axial stiffness K = (E A/L) b b^T, where b = [-n, n].
    class(truss_element_3d), intent(in) :: this
    !! The spatial truss element.
    integer(int32), intent(in), optional :: rule
    !! The integration rule is accepted for compatibility; stiffness is exact.
    real(real64), allocatable, dimension(:,:) :: rst
    !! The 6-by-6 global translational stiffness matrix.
    real(real64) :: length, direction(3), axial(6)
    integer(int32) :: row, col
    length = this%length()
    direction = [this%node_2%x - this%node_1%x, &
        this%node_2%y - this%node_1%y, this%node_2%z - this%node_1%z] / length
    axial = [-direction, direction]
    allocate(rst(6,6))
    do col = 1, 6
        do row = 1, 6
            rst(row,col) = this%material%modulus * this%area / length * axial(row) * axial(col)
        end do
    end do
end function

end module