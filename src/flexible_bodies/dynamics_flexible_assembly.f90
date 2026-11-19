module dynamics_flexible_assembly
    !! Connectivity-based assembly for flexible finite elements.
    !!
    !! Elements supply local operators in node-block order. The mesh gathers
    !! displacements and sums local contributions using a numbering derived
    !! from connectivity, independently of any stored element equation IDs.
    !! All outputs are unconstrained: apply boundary conditions separately.
    use iso_fortran_env, only : real64, int64
    use dynamics_flexible_body_types, only : fe_node
    use dynamics_finite_element, only : finite_element
    use linalg, only : csr_matrix, create_csr_matrix, create_empty_csr_matrix
    implicit none
    private

    public :: finite_element_mesh
    public :: mesh_element

    integer, parameter :: max_batch = 256
        !! Bound temporary element storage independently of mesh size.

    type :: mesh_element
        !! Owns one polymorphic element, allowing different types in one mesh.
        class(finite_element), allocatable :: element
            !! Allocate from a concrete element before assembly. Any pointer
            !! components of that element (such as a beam's material) still
            !! refer to their original targets, which must remain valid.
    end type

    type :: finite_element_mesh
        !! Owns a node table and a heterogeneous collection of elements.
        !! Shared nodes must use the same local DOF count and component order.
        !! The count is checked; compatible component order is the caller's
        !! responsibility. Numbering is recomputed on each call.
        type(fe_node), allocatable, dimension(:) :: nodes
            !! Reference nodes in equation-numbering order. Connectivity uses
            !! one-based positions, not node IDs. Remap element connectivity
            !! if the table is reordered.
        type(mesh_element), allocatable, dimension(:) :: elements
            !! Collection of concrete elements. An unallocated collection is
            !! empty; an allocated collection must have every slot populated.
    contains
        procedure :: number_dofs
            !! Return global DOF offsets by node, without evaluating elements.
        procedure, private :: mass_dense
        procedure, private :: mass_csr
        generic :: assemble_mass => mass_dense, mass_csr
            !! Assemble the global mass matrix in the format of the output
            !! argument: allocatable real rank two or `csr_matrix`.
        procedure, private :: stiffness_dense
        procedure, private :: stiffness_csr
        generic :: assemble_stiffness => stiffness_dense, stiffness_csr
            !! Assemble the global tangent at the supplied displacements in
            !! the format of the output argument: dense or CSR, not both.
            !! Optional `nthreads` enables parallel element evaluation.
        procedure :: assemble_internal_force
            !! Assemble element internal forces in global equation order.
    end type

    type :: element_map
        !! Index list for one element: node positions during layout, or
        !! global equations in local DOF order during assembly.
        integer, allocatable, dimension(:) :: indices
    end type

contains
! ------------------------------------------------------------------------------
    subroutine number_dofs(this, offsets)
        !! Return prefix sums of DOFs in node-table order. Node i owns
        !! equations offsets(i)+1 through offsets(i+1); unused nodes own none.
        !! Invalid connectivity or inconsistent DOF counts cause `error stop`.
        class(finite_element_mesh), intent(in) :: this
            !! Mesh whose topology determines the numbering.
        integer, allocatable, intent(out), dimension(:) :: offsets
            !! Allocated array of length `node_count + 1`, starting at zero
            !! and ending at the global DOF count. An empty mesh returns `[0]`.
        call build_layout(this, offsets)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine build_layout(this, offsets, maps)
        !! Validate topology, assign one contiguous block per referenced node,
        !! and optionally map each element's local DOFs to global equations.
        !! Reject missing elements, invalid or repeated node positions,
        !! inconsistent DOF counts, and unrepresentable global counts.
        class(finite_element_mesh), intent(in) :: this
        integer, allocatable, intent(out), dimension(:) :: offsets
            !! Prefix sums in node-table order; the final entry is the total.
        type(element_map), allocatable, intent(out), optional, dimension(:) :: maps
            !! One local-to-global equation map per element, when requested.

        integer, allocatable, dimension(:) :: connectivity, counts
        type(element_map), allocatable, dimension(:) :: topology
        integer :: nnode, nelem, i, j, k, node, nd, first

        nnode = 0
        nelem = 0
        if (allocated(this%nodes)) nnode = size(this%nodes)
        if (allocated(this%elements)) nelem = size(this%elements)
        allocate(counts(nnode), source=0)
        allocate(topology(nelem))

        do i = 1, nelem
            if (.not.allocated(this%elements(i)%element)) error stop "Mesh element is not allocated"
            associate (element => this%elements(i)%element)
                connectivity = element%node_connectivity()
                nd = element%dof_per_node()
                if (size(connectivity) == 0 .or. size(connectivity) /= element%node_count() .or. nd <= 0) &
                    error stop "Invalid mesh element topology"
                if (int(size(connectivity), int64)*nd /= element%dof_count()) &
                    error stop "Mesh element DOFs must be grouped by node"
                do j = 1, size(connectivity)
                    node = connectivity(j)
                    if (node < 1 .or. node > nnode) error stop "Mesh element node index out of bounds"
                    if (any(connectivity(:j-1) == node)) error stop "Mesh element repeats a node"
                    if (counts(node) /= 0 .and. counts(node) /= nd) &
                        error stop "Mesh elements disagree on the DOF count at a shared node"
                    counts(node) = nd
                end do
            end associate
            call move_alloc(connectivity, topology(i)%indices)
        end do

        allocate(offsets(nnode+1))
        offsets(1) = 0
        do i = 1, nnode
            if (counts(i) > huge(0) - 1 - offsets(i)) error stop "Mesh DOF count exceeds integer range"
            offsets(i+1) = offsets(i) + counts(i)
        end do
        if (.not.present(maps)) return

        allocate(maps(nelem))
        do i = 1, nelem
            nd = this%elements(i)%element%dof_per_node()
            allocate(maps(i)%indices(this%elements(i)%element%dof_count()))
            do j = 1, size(topology(i)%indices)
                node = topology(i)%indices(j)
                first = (j-1)*nd
                maps(i)%indices(first+1:first+nd) = [(offsets(node)+k, k=1,nd)]
            end do
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine mass_dense(this, matrix)
        !! Assemble mass into a newly allocated global dense matrix.
        class(finite_element_mesh), intent(inout) :: this
            !! Mesh; element routines may update their workspaces.
        real(real64), allocatable, intent(out), dimension(:,:) :: matrix
            !! Square matrix with one row and column per global DOF.
        call assemble_matrix(this, dense=matrix)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine mass_csr(this, matrix)
        !! Assemble mass into CSR without allocating a dense global matrix.
        class(finite_element_mesh), intent(inout) :: this
            !! Mesh; element routines may update their workspaces.
        type(csr_matrix), intent(out) :: matrix
            !! Square CSR matrix with sorted, unique column indices per row.
        call assemble_matrix(this, sparse=matrix)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine stiffness_dense(this, q, matrix, nthreads)
        !! Assemble tangent stiffness into a newly allocated dense matrix.
        class(finite_element_mesh), intent(inout) :: this
            !! Mesh; element routines may update their workspaces.
        real(real64), intent(in), dimension(:) :: q
            !! Displacements in global equation order. Wrong length is an error.
        real(real64), allocatable, intent(out), dimension(:,:) :: matrix
            !! Square matrix with one row and column per global DOF.
        integer, intent(in), optional :: nthreads
            !! Positive thread limit (default 1). Values above 1 require thread-safe element evaluation.
        call assemble_matrix(this, q=q, dense=matrix, nthreads=nthreads)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine stiffness_csr(this, q, matrix, nthreads)
        !! Assemble tangent stiffness into CSR without a dense global matrix.
        class(finite_element_mesh), intent(inout) :: this
            !! Mesh; element routines may update their workspaces.
        real(real64), intent(in), dimension(:) :: q
            !! Displacements in global equation order. Wrong length is an error.
        type(csr_matrix), intent(out) :: matrix
            !! Square CSR matrix with sorted, unique column indices per row.
            !! The sparsity pattern may change with `q`.
        integer, intent(in), optional :: nthreads
            !! Positive thread limit (default 1). Values above 1 require thread-safe element evaluation.
        call assemble_matrix(this, q=q, sparse=matrix, nthreads=nthreads)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine assemble_matrix(this, q, dense, sparse, nthreads)
        !! Assemble mass when `q` is absent, or tangent stiffness otherwise.
        !! Exactly one output is required. Dense assembly adds entries in
        !! place; sparse assembly collects nonzero triplets and reduces them
        !! to unique CSR entries after sorting.
        class(finite_element_mesh), intent(inout) :: this
        real(real64), intent(in), optional, dimension(:) :: q
            !! Global displacements for tangent assembly.
        real(real64), allocatable, intent(out), optional, dimension(:,:) :: dense
            !! Dense output in global equation order.
        type(csr_matrix), intent(out), optional :: sparse
            !! CSR output in the same equation order. Exact zeros are omitted,
            !! including entries that cancel when contributions are summed.
        integer, intent(in), optional :: nthreads

        type(element_map), allocatable, dimension(:) :: maps
        real(real64), allocatable, dimension(:,:,:) :: work
        integer, allocatable, dimension(:) :: rows, cols, offsets
        real(real64), allocatable, dimension(:) :: values
        integer(int64) :: capacity
        integer :: n, nd, max_nd, pos, i, j, k, threads, batch, first, last

        if (present(dense) .eqv. present(sparse)) error stop "Select exactly one mesh matrix format"
        call build_layout(this, offsets, maps)
        n = offsets(size(offsets))
        if (present(q)) then
            if (size(q) /= n) error stop "Mesh displacement vector has the wrong length"
        end if
        if (present(dense)) allocate(dense(n,n), source=0.0d0)
        max_nd = 0
        capacity = 0
        do i = 1, size(maps)
            nd = size(maps(i)%indices)
            max_nd = max(max_nd, nd)
            if (present(sparse)) then
                capacity = capacity + int(nd, int64)*nd
                if (capacity >= huge(pos)) error stop "Mesh sparse assembly exceeds integer range"
            end if
        end do
        if (present(sparse)) allocate(rows(capacity), cols(capacity), values(capacity))
        threads = assembly_threads(nthreads, size(maps))
        batch = 1
        if (threads > 1) batch = min(max_batch, size(maps))
        allocate(work(max_nd,max_nd,batch))

        pos = 0
        first = 1
        do while (first <= size(maps))
            last = first + min(batch - 1, size(maps) - first)
            if (threads > 1 .and. last > first) then
                !$omp parallel do default(shared) private(nd) schedule(static) num_threads(min(threads,last-first+1))
                do i = first, last
                    associate (indices => maps(i)%indices)
                        nd = size(indices)
                        call this%elements(i)%element%stiffness_matrix( &
                            q(indices), work(:nd,:nd,i-first+1), this%nodes)
                    end associate
                end do
                !$omp end parallel do
            else
                associate (indices => maps(first)%indices, element => this%elements(first)%element)
                    nd = size(indices)
                    if (present(q)) then
                        call element%stiffness_matrix(q(indices), work(:nd,:nd,1), this%nodes)
                    else
                        call element%mass_matrix(work(:nd,:nd,1), this%nodes)
                    end if
                end associate
            end if

            ! Shared equations are accumulated in element order.
            do i = first, last
                associate (indices => maps(i)%indices, local => work(:,:,i-first+1))
                    nd = size(indices)
                    if (present(dense)) then
                        do j = 1, nd
                            do k = 1, nd
                                dense(indices(k),indices(j)) = dense(indices(k),indices(j)) + local(k,j)
                            end do
                        end do
                    else
                        do j = 1, nd
                            do k = 1, nd
                                if (local(k,j) == 0.0d0) cycle
                                pos = pos + 1
                                rows(pos) = indices(k)
                                cols(pos) = indices(j)
                                values(pos) = local(k,j)
                            end do
                        end do
                    end if
                end associate
            end do
            first = last + 1
        end do
        if (present(sparse)) then
            if (pos == 0) then
                sparse = create_empty_csr_matrix(n, n, 0)
                sparse%row_indices = 1
            else
                sparse = create_csr_matrix(n, n, rows(:pos), cols(:pos), values(:pos))
                call sum_shared_entries(sparse)
            end if
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine sum_shared_entries(matrix)
        !! Coalesce consecutive duplicate columns in each sorted CSR row,
        !! drop exact cancellations, and update the row pointers in place.
        type(csr_matrix), intent(inout) :: matrix
            !! CSR matrix from the triplet converter, sorted but not coalesced.

        integer :: i, read_pos, write_pos, last, col
        real(real64) :: value

        read_pos = 1
        write_pos = 1
        do i = 1, size(matrix%row_indices)-1
            last = matrix%row_indices(i+1)-1
            matrix%row_indices(i) = write_pos
            do while (read_pos <= last)
                col = matrix%column_indices(read_pos)
                value = matrix%values(read_pos)
                read_pos = read_pos + 1
                do while (read_pos <= last)
                    if (matrix%column_indices(read_pos) /= col) exit
                    value = value + matrix%values(read_pos)
                    read_pos = read_pos + 1
                end do
                if (value == 0.0d0) cycle
                matrix%column_indices(write_pos) = col
                matrix%values(write_pos) = value
                write_pos = write_pos + 1
            end do
        end do
        matrix%row_indices(size(matrix%row_indices)) = write_pos
        matrix%column_indices = matrix%column_indices(:write_pos-1)
        matrix%values = matrix%values(:write_pos-1)
    end subroutine

! ------------------------------------------------------------------------------
    subroutine assemble_internal_force(this, q, f, nthreads)
        !! Gather local displacements and sum element forces at shared DOFs.
        !! This preserves the element force sign conventions and applies no
        !! external loads or boundary conditions.
        class(finite_element_mesh), intent(inout) :: this
            !! Mesh; element routines may update their workspaces.
        real(real64), intent(in), dimension(:) :: q
            !! Displacements in global equation order. Wrong length is an error.
        real(real64), allocatable, intent(out), dimension(:) :: f
            !! Allocated force vector with one entry per global DOF.
        integer, intent(in), optional :: nthreads
            !! Positive thread limit (default 1). Values above 1 require thread-safe element evaluation.

        type(element_map), allocatable, dimension(:) :: maps
        integer, allocatable, dimension(:) :: offsets
        real(real64), allocatable, dimension(:,:) :: work
        integer :: i, n, nd, max_nd, threads, batch, first, last

        call build_layout(this, offsets, maps)
        n = offsets(size(offsets))
        if (size(q) /= n) error stop "Mesh displacement vector has the wrong length"
        allocate(f(n), source=0.0d0)
        max_nd = 0
        do i = 1, size(maps)
            max_nd = max(max_nd, size(maps(i)%indices))
        end do
        threads = assembly_threads(nthreads, size(maps))
        batch = 1
        if (threads > 1) batch = min(max_batch, size(maps))
        allocate(work(max_nd,batch))
        first = 1
        do while (first <= size(maps))
            last = first + min(batch - 1, size(maps) - first)
            if (threads > 1 .and. last > first) then
                !$omp parallel do default(shared) private(nd) schedule(static) num_threads(min(threads,last-first+1))
                do i = first, last
                    associate (indices => maps(i)%indices)
                        nd = size(indices)
                        call this%elements(i)%element%internal_force_vector( &
                            q(indices), work(:nd,i-first+1), this%nodes)
                    end associate
                end do
                !$omp end parallel do
            else
                associate (indices => maps(first)%indices)
                    nd = size(indices)
                    call this%elements(first)%element%internal_force_vector(q(indices), work(:nd,1), this%nodes)
                end associate
            end if
            do i = first, last
                associate (indices => maps(i)%indices)
                    f(indices) = f(indices) + work(:size(indices),i-first+1)
                end associate
            end do
            first = last + 1
        end do
    end subroutine

! ------------------------------------------------------------------------------
    integer function assembly_threads(nthreads, nelem) result(threads)
        !! Without OpenMP the conditional assignment is omitted.
        integer, intent(in), optional :: nthreads
        integer, intent(in) :: nelem

        threads = 1
        if (present(nthreads)) then
            if (nthreads < 1) error stop "Mesh assembly thread count must be positive"
            !$ threads = max(1, min(nthreads, nelem, max_batch))
        end if
    end function

! ------------------------------------------------------------------------------
end module