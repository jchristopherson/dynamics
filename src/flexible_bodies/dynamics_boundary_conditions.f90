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
module dynamics_boundary_conditions
    !! Applies force and displacement boundary conditions to the unconstrained
    !! equations of a finite_element_mesh. A condition addresses a local DOF of
    !! a mesh node. Force conditions are added to the global force vector.
    !! Displacement conditions are eliminated from the assembled system and can
    !! be restored after the reduced system is solved.
    use iso_fortran_env, only : int32, real64
    use dynamics_error_handling, only : DYN_INDEX_OUT_OF_RANGE, &
        DYN_ARRAY_SIZE_ERROR, DYN_MATRIX_SIZE_ERROR
    use dynamics_flexible_assembly, only : finite_element_mesh
    use dynamics_structural, only : apply_boundary_conditions, &
        restore_constrained_values
    use linalg, only : csr_matrix, matmul, size
    implicit none
    private

    public :: nodal_condition
    public :: apply_force_boundary_conditions
    public :: apply_displacement_boundary_conditions
    public :: restore_displacement_boundary_conditions

    type :: nodal_condition
        !! Prescribes one local DOF of one node in a finite_element_mesh.
        integer(int32) :: node = 0
            !! The one-based position of the node in the mesh node table.
        integer(int32) :: component = 0
            !! The one-based local DOF within that node's block of equations.
        real(real64) :: value = 0.0d0
            !! The force (for force conditions) or displacement (for
            !! displacement conditions) applied to the DOF.
    end type

    interface apply_displacement_boundary_conditions
        module procedure :: apply_displacement_dense
        module procedure :: apply_displacement_csr
    end interface

contains
! ------------------------------------------------------------------------------
    subroutine apply_force_boundary_conditions(mesh, conditions, f)
        !! Adds force conditions to the global force vector. Repeated conditions
        !! on the same DOF accumulate.
        type(finite_element_mesh), intent(in) :: mesh
            !! The mesh that defines the global DOF numbering.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! The applied forces.
        real(real64), allocatable, intent(inout), dimension(:) :: f
            !! The global force vector. If unallocated, it is allocated as zeros
            !! with one entry per global DOF. Otherwise, its size must match.

        ! Local Variables
        integer(int32) :: ndof, i
        integer(int32), allocatable, dimension(:) :: dofs

        ! Process
        call condition_layout(mesh, conditions, ndof, dofs)
        if (.not.allocated(f)) then
            allocate(f(ndof), source = 0.0d0)
        else
            if (size(f) /= ndof) error stop DYN_ARRAY_SIZE_ERROR
        end if
        do i = 1, size(dofs)
            f(dofs(i)) = f(dofs(i)) + conditions(i)%value
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine apply_displacement_dense(mesh, conditions, k, f)
        !! Eliminates prescribed displacements from a dense global stiffness
        !! matrix and force vector. The reduced force is f - K*u_p, with the
        !! prescribed rows and columns removed.
        type(finite_element_mesh), intent(in) :: mesh
            !! The mesh that defines the global DOF numbering.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! The prescribed displacements. Duplicate DOFs are not allowed.
        real(real64), allocatable, intent(inout), dimension(:,:) :: k
            !! The global stiffness matrix. On output, it is reduced to the
            !! unconstrained DOFs.
        real(real64), allocatable, intent(inout), dimension(:) :: f
            !! The global force vector. On output, it is reduced to the
            !! unconstrained DOFs.

        ! Local Variables
        integer(int32) :: ndof
        integer(int32), allocatable, dimension(:) :: dofs, gdof
        real(real64), allocatable, dimension(:) :: u

        ! Input Checking
        call condition_layout(mesh, conditions, ndof, dofs)
        if (.not.allocated(k)) error stop DYN_MATRIX_SIZE_ERROR
        if (size(k, 1) /= ndof .or. size(k, 2) /= ndof) error stop DYN_MATRIX_SIZE_ERROR
        if (.not.allocated(f)) error stop DYN_ARRAY_SIZE_ERROR
        if (size(f) /= ndof) error stop DYN_ARRAY_SIZE_ERROR

        ! Build the prescribed displacement vector and move its contribution to the force
        allocate(u(ndof), source = 0.0d0)
        u(dofs) = conditions%value
        f = f - matmul(k, u)

        ! Reduce both the matrix and the vector
        gdof = dofs
        k = apply_boundary_conditions(gdof, k)
        gdof = dofs
        f = apply_boundary_conditions(gdof, f)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine apply_displacement_csr(mesh, conditions, k, f)
        !! Eliminates prescribed displacements from a CSR global stiffness
        !! matrix and a force vector. The reduced force is f - K*u_p, with the
        !! prescribed rows and columns removed.
        type(finite_element_mesh), intent(in) :: mesh
            !! The mesh that defines the global DOF numbering.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! The prescribed displacements. Duplicate DOFs are not allowed.
        type(csr_matrix), intent(inout) :: k
            !! The global stiffness matrix. On output, it is reduced to the
            !! unconstrained DOFs.
        real(real64), allocatable, intent(inout), dimension(:) :: f
            !! The global force vector. On output, it is reduced to the
            !! unconstrained DOFs.

        ! Local Variables
        integer(int32) :: ndof
        integer(int32), allocatable, dimension(:) :: dofs, gdof
        real(real64), allocatable, dimension(:) :: u

        ! Input Checking
        call condition_layout(mesh, conditions, ndof, dofs)
        if (size(k, 1) /= ndof .or. size(k, 2) /= ndof) error stop DYN_MATRIX_SIZE_ERROR
        if (.not.allocated(f)) error stop DYN_ARRAY_SIZE_ERROR
        if (size(f) /= ndof) error stop DYN_ARRAY_SIZE_ERROR

        ! Build the prescribed displacement vector and move its contribution to the force
        allocate(u(ndof), source = 0.0d0)
        u(dofs) = conditions%value
        f = f - matmul(k, u)

        ! Reduce both the matrix and the vector
        gdof = dofs
        k = apply_boundary_conditions(gdof, k)
        gdof = dofs
        f = apply_boundary_conditions(gdof, f)
    end subroutine

! ------------------------------------------------------------------------------
    function restore_displacement_boundary_conditions(mesh, conditions, uf) result(u)
        !! Expands the solution of the reduced system to all global DOFs by
        !! reinserting the prescribed displacements.
        type(finite_element_mesh), intent(in) :: mesh
            !! The mesh that defines the global DOF numbering.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! The prescribed displacements used to reduce the system.
        real(real64), intent(in), dimension(:) :: uf
            !! The solution for the unconstrained DOFs.
        real(real64), allocatable, dimension(:) :: u
            !! The displacement vector for all global DOFs.

        ! Local Variables
        integer(int32) :: ndof
        integer(int32), allocatable, dimension(:) :: dofs, gdof

        ! Input Checking
        call condition_layout(mesh, conditions, ndof, dofs)
        if (size(uf) /= ndof - size(dofs)) error stop DYN_ARRAY_SIZE_ERROR

        ! Process
        gdof = dofs
        u = restore_constrained_values(gdof, uf)
        u(dofs) = conditions%value
    end function

! ------------------------------------------------------------------------------
    subroutine condition_layout(mesh, conditions, ndof, dofs)
        !! Maps each condition to its global DOF using the mesh numbering.
        type(finite_element_mesh), intent(in) :: mesh
            !! The mesh that defines the global DOF numbering.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! The conditions to map.
        integer(int32), intent(out) :: ndof
            !! The total number of global DOFs in the mesh.
        integer(int32), allocatable, intent(out), dimension(:) :: dofs
            !! The global DOF of each condition, in the order of conditions.

        ! Local Variables
        integer, allocatable, dimension(:) :: offsets
        integer(int32) :: i, n, c, nnode

        ! Initialization
        call mesh%number_dofs(offsets)
        ndof = int(offsets(size(offsets)), int32)
        nnode = int(size(offsets) - 1, int32)
        allocate(dofs(size(conditions)))

        ! Map each condition, rejecting nodes or components outside the mesh
        do i = 1, size(conditions)
            n = conditions(i)%node
            c = conditions(i)%component
            if (n < 1 .or. n > nnode) error stop DYN_INDEX_OUT_OF_RANGE
            if (c < 1 .or. c > offsets(n+1) - offsets(n)) error stop DYN_INDEX_OUT_OF_RANGE
            dofs(i) = offsets(n) + c
        end do
    end subroutine
end module