module dynamics_structural
    use iso_fortran_env
    use linalg, only : csr_matrix, create_csr_matrix, dense_to_csr, sort, &
        size, assignment(=), lu_factor, solve_lu, sparse_direct_solve
    use dynamics_error_handling
    use dynamics_geometry
    implicit none
    private
    public :: csr_matrix
    public :: assignment(=)
    public :: DYN_ONE_POINT_INTEGRATION_RULE
    public :: DYN_TWO_POINT_INTEGRATION_RULE
    public :: DYN_THREE_POINT_INTEGRATION_RULE
    public :: DYN_FOUR_POINT_INTEGRATION_RULE
    public :: node
    public :: material
    public :: element
    public :: shape_function_derivative
    public :: shape_function_second_derivative
    public :: create_connectivity_matrix
    public :: assemble_static_system
    public :: assemble_dynamic_system
    public :: apply_boundary_conditions
    public :: apply_displacement_constraint
    public :: restore_constrained_values
    public :: solve_static_system
    public :: line_element

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: DYN_ONE_POINT_INTEGRATION_RULE = 1
        !! Defines a single-point integration rule.
    integer(int32), parameter :: DYN_TWO_POINT_INTEGRATION_RULE = 2
        !! Defines a two-point integration rule.
    integer(int32), parameter :: DYN_THREE_POINT_INTEGRATION_RULE = 3
        !! Defines a three-point integration rule.
    integer(int32), parameter :: DYN_FOUR_POINT_INTEGRATION_RULE = 4
        !! Defines a four-point integration rule.

! ******************************************************************************
! TYPES
! ------------------------------------------------------------------------------
    type, extends(point) :: node
        !! Defines a node.
        integer(int32) :: index
            !! The global index of the node.
        integer(int32) :: dof
            !! The number of degrees of freeedom associated with this node.
    end type

    interface node
        module procedure :: nd_init_1
        module procedure :: nd_init_2
    end interface

! ------------------------------------------------------------------------------
    type :: material
        !! Defines a linear-elastic-isotropic material.
        real(real64) :: density
            !! The density of the material.
        real(real64) :: modulus
            !! The modulus of elasticity of the material.
        real(real64) :: poissons_ratio
            !! The Poisson's ratio of the material.
    end type

    interface material
        module procedure :: mat_init
    end interface

! ------------------------------------------------------------------------------
    type, abstract :: element
        !! Defines an element.
        type(material) :: material
            !! The material.
    contains
        procedure(element_query), deferred, public, pass :: get_dimensionality
        procedure(element_query), deferred, public, pass :: get_node_count
        procedure(element_get_node), deferred, public, pass :: get_node
        procedure(element_query), deferred, public, pass :: get_dof_per_node
        procedure(element_shape_function), deferred, public, pass :: &
            evaluate_shape_function
        procedure(element_matrix_function), deferred, public, pass :: &
            shape_function_matrix
        procedure(element_matrix_function), deferred, public, pass :: &
            strain_displacement_matrix
        procedure(element_const_matrix_function), deferred, public, &
            pass :: constitutive_matrix
        procedure(element_matrix_function), deferred, public, pass :: &
            jacobian
        procedure, public :: stiffness_matrix => e_stiffness_matrix
        procedure, public :: mass_matrix => e_mass_matrix
        procedure, public :: external_force_vector => e_ext_force_vector
    end type
    
! ------------------------------------------------------------------------------
    type, extends(element), abstract :: line_element
        !! Defines a line element type.
        real(real64) :: area
            !! The element cross-sectional area.
    contains
        procedure(line_element_get_terminal), deferred, public, pass :: &
            get_terminal_nodes
        procedure(line_element_const_matrix_function), deferred, public, &
            pass :: rotation_matrix
        procedure, public :: length => le_length
        procedure, public :: stiffness_matrix => le_stiffness_matrix
        procedure, public :: mass_matrix => le_mass_matrix
        procedure, public :: external_force_vector => le_ext_force_vector
    end type

! ******************************************************************************
! INTERFACES
! ------------------------------------------------------------------------------
    interface
        pure function element_query(this) result(rst)
            !! Defines the signature of a function performing a query on an
            !! integer-valued property of a element type.
            use iso_fortran_env, only : int32
            import element
            class(element), intent(in) :: this
                !! The element object.
            integer(int32) :: rst
                !! The resulting value.
        end function

        pure function element_get_node(this, i) result(rst)
            !! Defines the signature of a function for retrieving the requested
            !! node from the element.
            use iso_fortran_env, only : int32
            import node
            import element
            class(element), intent(in) :: this
                !! The element object.
            integer(int32), intent(in) :: i
                !! The local index of the node to retrieve.
            type(node) :: rst
                !! The node.
        end function

        pure function element_matrix_function(this, s) result(rst)
            !! Defines the signature of a routine for returning a matrix
            !! associated with the element.
            use iso_fortran_env, only : real64
            import element
            class(element), intent(in) :: this
                !! The element object.
            real(real64), intent(in), dimension(:) :: s
                !! The value of the natural coordinates at which the matrix
                !! should be evaluated.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The resulting matrix.
        end function

        pure function element_const_matrix_function(this) result(rst)
            !! Defines the signature of a routine for returning a matrix
            !! associated with the element.
            use iso_fortran_env, only : real64
            import element
            class(element), intent(in) :: this
                !! The element object.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The resulting matrix.
        end function

        pure function element_shape_function(this, i, s) result(rst)
            !! Defines the signature of a routine for computing the value of
            !! the i-th element shape function at natural coordinate.
            use iso_fortran_env, only : int32, real64
            import element
            class(element), intent(in) :: this
                !! The element object.
            integer(int32), intent(in) :: i
                !! The index of the shape function to evaluate.
            real(real64), intent(in), dimension(:) :: s
                !! The value of the natural coordinates at which to evaluate
                !! the shape function.
            real(real64) :: rst
                !! The value of the i-th shape function at s.
        end function
        
        pure function integrand(elem, s) result(rst)
            !! Defines the signature of a function containing an integrand.
            use iso_fortran_env, only : real64
            import element
            class(element), intent(in) :: elem
                !! The element object.
            real(real64), intent(in), dimension(:) :: s
                !! The natural coordinate at which to evaluate the integrand.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The result.
        end function

        pure subroutine line_element_get_terminal(this, i1, i2)
            !! Defines the signature of a routine for returning the terminal
            !! node numbers.
            use iso_fortran_env, only : int32
            import line_element
            class(line_element), intent(in) :: this
                !! The line_element object.
            integer(int32), intent(out) :: i1
                !! The index of the node at the head of the element.
            integer(int32), intent(out) :: i2
                !! The index of the node at the tail of the element.
        end subroutine

        pure function line_element_const_matrix_function(this) result(rst)
            !! Defines the signature of a routine for returning a matrix
            !! associated with the line_element.
            use iso_fortran_env, only : real64
            import line_element
            class(line_element), intent(in) :: this
                !! The line_element object.
            real(real64), allocatable, dimension(:,:) :: rst
                !! The resulting matrix.
        end function
    end interface

! ******************************************************************************
! OVERLOADED ROUTINES
! ------------------------------------------------------------------------------
    interface apply_boundary_conditions
        module procedure :: apply_boundary_conditions_mtx
        module procedure :: apply_boundary_conditions_vec
        module procedure :: apply_boundary_conditions_csr
    end interface

    interface restore_constrained_values
        module procedure :: restore_constrained_values_dense
        module procedure :: restore_constrained_values_csr
    end interface

    interface apply_displacement_constraint
        module procedure :: apply_displacement_constraint_dense
        module procedure :: apply_displacement_constraint_csr
    end interface

    interface assemble_static_system
        module procedure :: assemble_static_system_dense
        module procedure :: assemble_static_system_csr
    end interface

    interface assemble_dynamic_system
        module procedure :: assemble_dynamic_system_dense
        module procedure :: assemble_dynamic_system_csr
    end interface

    interface solve_static_system
        module procedure :: solve_static_system_dense
        module procedure :: solve_static_system_csr
    end interface
contains
! ******************************************************************************
! DIFFERENTIATION ROUTINES
! ------------------------------------------------------------------------------
pure function shape_function_derivative(index, elem, s, i) result(rst)
    !! Computes the derivative of the shape function with respect to the natural
    !! coordinate specified.
    !! The derivative is approximated centrally as
    !! $$ N_{,i}(\boldsymbol{s})\approx
    !! \frac{N(\boldsymbol{s}+h\boldsymbol{e}_i)-
    !! N(\boldsymbol{s}-h\boldsymbol{e}_i)}{2h}. $$
    !! The second derivative uses the centered finite difference
    !! $$ N_{,ii}(\boldsymbol{s})\approx
    !! \frac{N(\boldsymbol{s}+h\boldsymbol{e}_i)-2N(\boldsymbol{s})+
    !! N(\boldsymbol{s}-h\boldsymbol{e}_i)}{h^2}. $$
    integer(int32), intent(in) :: index
        !! The index of the shape function to evaluate.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the derivative.
    integer(int32), intent(in) :: i
        !! The index of the natural coordinate to with which the derivative is
        !! to be computed.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    real(real64) :: na, nb, h(size(s))

    ! Initialization
    h = 0.0d0
    h(i) = sqrt(epsilon(na))

    ! Process
    na = elem%evaluate_shape_function(index, s + h)
    nb = elem%evaluate_shape_function(index, s - h)
    rst = (na - nb) / (2.0d0 * h(i))
end function

! ------------------------------------------------------------------------------
pure function shape_function_second_derivative(index, elem, s, i) result(rst)
    !! Computes the second derivative of the shape function with respect to the
    !! natural coordinate specified.
    integer(int32), intent(in) :: index
        !! The index of the shape function to evaluate.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate at which to evaluate the derivative.
    integer(int32), intent(in) :: i
        !! The index of the natural coordinate to with which the derivative is
        !! to be computed.
    real(real64) :: rst
        !! The result.

    ! Local Variables
    real(real64) :: na, nb, nc, h(size(s))

    ! Initialization
    h = 0.0d0
    h(i) = (epsilon(na))**0.25d0
    na = elem%evaluate_shape_function(index, s + h)
    nb = elem%evaluate_shape_function(index, s)
    nc = elem%evaluate_shape_function(index, s - h)
    rst = (na - 2.0d0 * nb + nc) / (h(i)**2)
end function

! ******************************************************************************
! INTEGRATION
! ------------------------------------------------------------------------------
pure function get_model_parameters(rule) result(rst)
    !! Gets the requested integration model parameters.
    !! Each returned row contains a Gauss point and weight \((s_i,w_i)\) for
    !! approximating
    !! $$ \int_{-1}^{1}f(s)\,ds\approx\sum_i w_i f(s_i). $$
    integer(int32), intent(in) :: rule
        !! The integration rule.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integration parameters.

    ! Local Variables
    real(real64) :: x, w1, w2, pt1, pt2

    ! Process
    select case (rule)
    case (DYN_ONE_POINT_INTEGRATION_RULE)
        allocate(rst(2, 2))
        rst = reshape([0.0d0, 0.0d0, 2.0d0, 2.0d0], [2, 2])
    case (DYN_TWO_POINT_INTEGRATION_RULE)
        allocate(rst(2, 2))
        x = sqrt(3.0d0) / 3.0d0
        rst = reshape([-x, x, 1.0d0, 1.0d0], [2, 2])
    case (DYN_THREE_POINT_INTEGRATION_RULE)
        allocate(rst(3, 2))
        x = sqrt(3.0d0 / 5.0d0)
        rst = reshape([0.0d0, -x, x, w1, w2, w2], [3, 2])
    case default ! Four Point Rule
        allocate(rst(4, 2))
        pt1 = sqrt((3.0d0 / 7.0d0) - (2.0d0 / 7.0d0) * sqrt(6.0d0 / 5.0d0))
        pt1 = sqrt((3.0d0 / 7.0d0) + (2.0d0 / 7.0d0) * sqrt(6.0d0 / 5.0d0))
        w1 = (1.8d1 + sqrt(3.0d1)) / 3.6d1
        w2 = (1.8d1 - sqrt(3.0d1)) / 3.6d1
        rst = reshape([-pt1, pt1, -pt2, pt2, w1, w1, w2, w2], [4, 2])
    end select
end function

! ------------------------------------------------------------------------------
pure function integrate_1d(fcn, elem, rule) result(rst)
    !! Computes the integral of the specified integrand given an element and an
    !! integration rule.
    !! The element integral is evaluated by the Gauss rule in the element's
    !! natural coordinate \(s\in[-1,1]\).
    procedure(integrand) :: fcn
        !! The integrand.
    class(element), intent(in) :: elem
        !! The element object.
    integer(int32), intent(in) :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
    real(real64), allocatable, dimension(:,:) :: rst
        !! The result of the integration.

    ! Local Variables
    integer(int32) :: i
    real(real64), allocatable, dimension(:,:) :: s

    ! Process
    s = get_model_parameters(rule)
    rst = s(1,2) * fcn(elem, [s(1,1)])
    do i = 2, size(s, 1)
        rst = rst + s(i,2) * fcn(elem, [s(i,1)])
    end do
end function

! ------------------------------------------------------------------------------
pure function integrate(fcn, elem, rule) result(rst)
    !! Computes the integral of the specified integrand given an element and an
    !! integration rule.
    !! For an isoparametric element, the physical-coordinate integral includes
    !! the Jacobian determinant,
    !! \(\int_{\Omega_e}g\,d\Omega=\int_{-1}^{1}g(s)J(s)\,ds\).
    procedure(integrand) :: fcn
        !! The integrand.
    class(element), intent(in) :: elem
        !! The element object.
    integer(int32), intent(in) :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
    real(real64), allocatable, dimension(:,:) :: rst
        !! The result of the integration.

    ! Process
    select type (elem)
    class is (line_element)
        rst = integrate_1d(fcn, elem, rule)
    end select
end function

! ******************************************************************************
! ASSEMBLY ROUTINES
! ------------------------------------------------------------------------------
pure function find_global_dof(n, nodes) result(rst)
    !! Finds the index of the global DOF node in a list of nodes.
    class(node), intent(in) :: n
        !! The node for which to search.
    class(node), intent(in), dimension(:) :: nodes
        !! The list of nodes
    integer(int32) :: rst
        !! The requested index.

    ! Local Variables
    integer(int32) :: i

    ! Process
    rst = 0
    do i = 1, size(nodes)
        if (n%index == nodes(i)%index) then
            rst = rst + 1
            exit
        end if
        rst = rst + nodes(i)%dof
    end do
end function

! ------------------------------------------------------------------------------
function create_connectivity_matrix(gdof, e, nodes) result(rst)
    !! Creates a connectivity matrix for the element, stored in CSR format.
    !! The matrix contains exactly one non-zero (unity) entry per row;
    !! therefore, it is well-suited to a sparse representation.
    integer(int32), intent(in) :: gdof
        !! The number of global degrees of freedom.
    class(element), intent(in) :: e
        !! The element.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    type(csr_matrix) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: i, j, col, nnodes, nnz, row
    integer(int32), allocatable, dimension(:) :: rows, cols
    real(real64), allocatable, dimension(:) :: vals
    
    ! Initialization
    nnodes = e%get_node_count()
    nnz = e%get_dof_per_node() * nnodes
    allocate(rows(nnz), cols(nnz))
    allocate(vals(nnz), source = 1.0d0)

    ! Process
    row = 0
    do j = 1, nnodes
        col = find_global_dof(e%get_node(j), nodes)
        do i = 1, e%get_dof_per_node()
            row = row + 1
            rows(row) = row
            cols(row) = col
            col = col + 1
        end do
    end do
    rst = create_csr_matrix(nnz, gdof, rows, cols, vals)
end function

! ------------------------------------------------------------------------------
subroutine assemble_static_system_csr(gdof, elements, nodes, k, rule)
    !! Assembles the global stiffness matrix in CSR format.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(element), intent(in) :: elements(:)
        !! The finite elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    type(csr_matrix), intent(out) :: k
        !! The assembled global stiffness matrix in CSR format.
    integer(int32), intent(in), optional :: rule
        !! The numerical integration rule.

    ! Local Variables
    integer(int32) :: i, j, eidx, row, col, ndof
    real(real64), allocatable :: kdense(:,:), ke(:,:)

    ! Initialization
    allocate(kdense(gdof,gdof), source = 0.0d0)

    ! Accumulate element contributions in global work storage.
    do eidx = 1, size(elements)
        if (present(rule)) then
            ke = elements(eidx)%stiffness_matrix(rule)
        else
            ke = elements(eidx)%stiffness_matrix()
        end if
        ndof = size(ke, 1)
        do i = 1, ndof
            row = find_global_dof(elements(eidx)%get_node( &
                (i - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                mod(i - 1, elements(eidx)%get_dof_per_node())
            do j = 1, ndof
                col = find_global_dof(elements(eidx)%get_node( &
                    (j - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                    mod(j - 1, elements(eidx)%get_dof_per_node())
                kdense(row, col) = kdense(row, col) + ke(i, j)
            end do
        end do
    end do
    k = dense_to_csr(kdense)
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_dynamic_system_csr(gdof, elements, nodes, m, k, rule)
    !! Assembles global mass and stiffness matrices in CSR format.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(element), intent(in) :: elements(:)
        !! The finite elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    type(csr_matrix), intent(out) :: m
        !! The assembled global mass matrix in CSR format.
    type(csr_matrix), intent(out) :: k
        !! The assembled global stiffness matrix in CSR format.
    integer(int32), intent(in), optional :: rule
        !! The numerical integration rule.

    ! Local Variables
    integer(int32) :: i, j, eidx, row, col, ndof
    real(real64), allocatable :: mdense(:,:), kdense(:,:), km(:,:), ke(:,:)

    ! Initialization
    allocate(mdense(gdof,gdof), kdense(gdof,gdof), source = 0.0d0)

    ! Accumulate element contributions in global work storage.
    do eidx = 1, size(elements)
        if (present(rule)) then
            km = elements(eidx)%mass_matrix(rule)
            ke = elements(eidx)%stiffness_matrix(rule)
        else
            km = elements(eidx)%mass_matrix()
            ke = elements(eidx)%stiffness_matrix()
        end if
        ndof = size(ke, 1)
        do i = 1, ndof
            row = find_global_dof(elements(eidx)%get_node( &
                (i - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                mod(i - 1, elements(eidx)%get_dof_per_node())
            do j = 1, ndof
                col = find_global_dof(elements(eidx)%get_node( &
                    (j - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                    mod(j - 1, elements(eidx)%get_dof_per_node())
                mdense(row, col) = mdense(row, col) + km(i, j)
            end do
        end do
        do i = 1, ndof
            row = find_global_dof(elements(eidx)%get_node( &
                (i - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                mod(i - 1, elements(eidx)%get_dof_per_node())
            do j = 1, ndof
                col = find_global_dof(elements(eidx)%get_node( &
                    (j - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                    mod(j - 1, elements(eidx)%get_dof_per_node())
                kdense(row, col) = kdense(row, col) + ke(i, j)
            end do
        end do
    end do
    m = dense_to_csr(mdense)
    k = dense_to_csr(kdense)
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_static_system_dense(gdof, elements, nodes, k, rule)
    !! Assembles a dense global stiffness matrix.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(element), intent(in) :: elements(:)
        !! The finite elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    real(real64), allocatable, intent(out) :: k(:,:)
        !! The assembled global stiffness matrix.
    integer(int32), intent(in), optional :: rule
        !! The numerical integration rule.

    ! Local Variables
    integer(int32) :: i, j, eidx, row, col, ndof
    real(real64), allocatable :: ke(:,:)

    ! Initialization
    allocate(k(gdof, gdof), source = 0.0d0)

    ! Accumulate element contributions in global dense storage.
    do eidx = 1, size(elements)
        if (present(rule)) then
            ke = elements(eidx)%stiffness_matrix(rule)
        else
            ke = elements(eidx)%stiffness_matrix()
        end if
        ndof = size(ke, 1)
        do i = 1, ndof
            row = find_global_dof(elements(eidx)%get_node( &
                (i - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                mod(i - 1, elements(eidx)%get_dof_per_node())
            do j = 1, ndof
                col = find_global_dof(elements(eidx)%get_node( &
                    (j - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                    mod(j - 1, elements(eidx)%get_dof_per_node())
                k(row, col) = k(row, col) + ke(i, j)
            end do
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine assemble_dynamic_system_dense(gdof, elements, nodes, m, k, rule)
    !! Assembles dense global mass and stiffness matrices.
    integer(int32), intent(in) :: gdof
        !! The total number of global degrees of freedom.
    class(element), intent(in) :: elements(:)
        !! The finite elements to assemble.
    class(node), intent(in), dimension(:) :: nodes
        !! The global node list.
    real(real64), allocatable, intent(out) :: m(:,:)
        !! The assembled global mass matrix.
    real(real64), allocatable, intent(out) :: k(:,:)
        !! The assembled global stiffness matrix.
    integer(int32), intent(in), optional :: rule
        !! The numerical integration rule.

    ! Local Variables
    integer(int32) :: i, j, eidx, row, col, ndof
    real(real64), allocatable :: km(:,:), ke(:,:)

    ! Initialization
    allocate(m(gdof, gdof), k(gdof, gdof), source = 0.0d0)

    ! Accumulate element contributions in global dense storage.
    do eidx = 1, size(elements)
        if (present(rule)) then
            km = elements(eidx)%mass_matrix(rule)
            ke = elements(eidx)%stiffness_matrix(rule)
        else
            km = elements(eidx)%mass_matrix()
            ke = elements(eidx)%stiffness_matrix()
        end if
        ndof = size(ke, 1)
        do i = 1, ndof
            row = find_global_dof(elements(eidx)%get_node( &
                (i - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                mod(i - 1, elements(eidx)%get_dof_per_node())
            do j = 1, ndof
                col = find_global_dof(elements(eidx)%get_node( &
                    (j - 1) / elements(eidx)%get_dof_per_node() + 1), nodes) + &
                    mod(j - 1, elements(eidx)%get_dof_per_node())
                m(row, col) = m(row, col) + km(i, j)
                k(row, col) = k(row, col) + ke(i, j)
            end do
        end do
    end do
end subroutine

! ******************************************************************************
! BOUNDARY CONDITIONS ROUTINES
! ------------------------------------------------------------------------------
function apply_boundary_conditions_mtx(gdof, x) result(rst)
    !! Applies boundary conditions to a matrix by removal of the appropriate
    !! rows and columns.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    real(real64), intent(in), dimension(:,:) :: x
        !! The matrix to constrain.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The altered matrix.

    ! Local Variables
    integer(int32) :: i, j, ii, m, n, nbc, mnew
    integer(int32), allocatable, dimension(:) :: indices
    
    ! Initialization
    m = size(x, 1)
    n = size(x, 2)
    nbc = size(gdof)
    mnew = m - nbc

    ! Input Checking
    if (m /= n) error stop DYN_MATRIX_SIZE_ERROR
    if (mnew < 1) error stop DYN_CONSTRAINT_ERROR
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > m) error stop DYN_INDEX_OUT_OF_RANGE
    end do

    ! Memory Allocation
    allocate(rst(mnew, mnew), indices(m - nbc))

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) error stop DYN_NONMONOTONIC_ARRAY_ERROR
    end do

    ! Process
    ii = 1
    j = 0
    do i = 1, m
        if (gdof(ii) /= i) then
            j = j + 1
            indices(j) = i
        else
            ii = ii + 1
            if (ii > nbc) ii = nbc
        end if
    end do

    ! Now, we only need store the rows and columns stored in indices
    rst = x(indices,indices)
end function

! ------------------------------------------------------------------------------
function apply_boundary_conditions_csr(gdof, x) result(rst)
    !! Applies boundary conditions to a CSR-format sparse matrix by removal of
    !! the appropriate rows and columns.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    type(csr_matrix), intent(in) :: x
        !! The matrix to constrain.
    type(csr_matrix) :: rst
        !! The altered matrix.

    ! Local Variables
    integer(int32) :: i, ii, j, k, m, n, nbc, mnew, nnz
    integer(int32), allocatable, dimension(:) :: map, rows, cols
    real(real64), allocatable, dimension(:) :: vals
    
    ! Initialization
    m = size(x, 1)
    n = size(x, 2)
    nbc = size(gdof)
    mnew = m - nbc

    ! Input Checking
    if (m /= n) error stop DYN_MATRIX_SIZE_ERROR
    if (mnew < 1) error stop DYN_CONSTRAINT_ERROR
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > m) error stop DYN_INDEX_OUT_OF_RANGE
    end do

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) error stop DYN_NONMONOTONIC_ARRAY_ERROR
    end do

    ! Build a map from the old row/column index to the new, constrained index;
    ! a value of zero denotes a row/column that is to be removed
    allocate(map(m))
    ii = 1
    j = 0
    do i = 1, m
        if (ii <= nbc) then
            if (gdof(ii) == i) then
                map(i) = 0
                ii = ii + 1
                cycle
            end if
        end if
        j = j + 1
        map(i) = j
    end do

    ! Count the number of retained non-zero entries
    nnz = 0
    do i = 1, m
        if (map(i) == 0) cycle
        do k = x%row_indices(i), x%row_indices(i+1) - 1
            if (map(x%column_indices(k)) == 0) cycle
            nnz = nnz + 1
        end do
    end do

    ! Populate the retained entries, remapped to the new index set
    allocate(rows(nnz), cols(nnz), vals(nnz))
    nnz = 0
    do i = 1, m
        if (map(i) == 0) cycle
        do k = x%row_indices(i), x%row_indices(i+1) - 1
            j = x%column_indices(k)
            if (map(j) == 0) cycle
            nnz = nnz + 1
            rows(nnz) = map(i)
            cols(nnz) = map(j)
            vals(nnz) = x%values(k)
        end do
    end do
    rst = create_csr_matrix(mnew, mnew, rows, cols, vals)
end function

! ------------------------------------------------------------------------------
function apply_boundary_conditions_vec(gdof, x) result(rst)
    !! Applies boundary conditions to a vector by removal of the appropriate
    !! items.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    real(real64), intent(in), dimension(:) :: x
        !! The vector to constrain.
    real(real64), allocatable, dimension(:) :: rst
        !! The altered vector.

    ! Local Variables
    integer(int32) :: i, j, ii, n, nbc, nnew
    integer(int32), allocatable, dimension(:) :: indices
    
    ! Initialization
    n = size(x)
    nbc = size(gdof)
    nnew = n - nbc

    ! Input Checking
    if (nnew < 1) error stop DYN_CONSTRAINT_ERROR
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > n) error stop DYN_INDEX_OUT_OF_RANGE
    end do

    ! Memory Allocation
    allocate(rst(nnew), indices(n - nbc))

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) error stop DYN_NONMONOTONIC_ARRAY_ERROR
    end do

    ! Process
    ii = 1
    j = 0
    do i = 1, n
        if (gdof(ii) /= i) then
            j = j + 1
            indices(j) = i
        else
            ii = ii + 1
            if (ii > nbc) ii = nbc
        end if
    end do

    ! Now, just store the appropriate items in the output vector
    rst = x(indices)
end function

! ------------------------------------------------------------------------------
function restore_constrained_values_dense(gdof, x) result(rst)
    !! Restores the constrained degrees-of-freedom from the boundary conditions
    !! applied by apply_boundary_conditions.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    real(real64), intent(in), dimension(:) :: x
        !! The constrained vector.
    real(real64), allocatable, dimension(:) :: rst
        !! The altered vector.

    ! Local Variables
    integer(int32) ::i, j, ii, n, nbc, nnew
    
    ! Initialization
    n = size(x)
    nbc = size(gdof)
    nnew = n + nbc

    ! Input Checking
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > nnew) error stop DYN_INDEX_OUT_OF_RANGE
    end do

    ! Memory Allocation
    allocate(rst(nnew), source = 0.0d0)

    ! Sort gdof into ascending order
    call sort(gdof, .true.)

    ! Check for duplicate values in GDOF
    do i = 2, nbc
        if (gdof(i) == gdof(i-1)) error stop DYN_NONMONOTONIC_ARRAY_ERROR
    end do

    ! Process
    ii = 1
    j = 0
    do i = 1, nnew
        if (i == gdof(ii)) then
            ii = ii + 1
            if (ii > nbc) ii = nbc
        else
            j = j + 1
            rst(i) = x(j)
        end if
    end do
end function

! ------------------------------------------------------------------------------
function restore_constrained_values_csr(gdof, x) result(rst)
    !! Restores constrained rows and columns to a reduced CSR matrix.
    integer(int32), intent(inout), dimension(:) :: gdof
        !! An array of the global degrees of freedom to restrain.  The array
        !! is sorted into ascending order on output.
    type(csr_matrix), intent(in) :: x
        !! The reduced CSR matrix.
    type(csr_matrix) :: rst
        !! The expanded CSR matrix with zero constrained rows and columns.

    ! Local Variables
    integer(int32) :: i, j, ii, n, nbc, nnew, nnz, pos, nout
    integer(int32), allocatable :: indices(:), rows(:), cols(:)
    real(real64), allocatable :: vals(:)

    ! Initialization
    n = size(x, 1)
    nbc = size(gdof)
    nnew = n + nbc
    nnz = size(x%values)

    ! Input Checking
    if (size(x, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
    do i = 1, nbc
        if (gdof(i) < 1 .or. gdof(i) > nnew) error stop DYN_INDEX_OUT_OF_RANGE
    end do

    ! Build the map from reduced indices to unconstrained global indices.
    allocate(indices(n))
    nout = 0
    if (nbc == 0) then
        do i = 1, n
            nout = nout + 1
            indices(nout) = i
        end do
    else
        call sort(gdof, .true.)
        do i = 2, nbc
            if (gdof(i) == gdof(i-1)) &
                error stop DYN_NONMONOTONIC_ARRAY_ERROR
        end do
        ii = 1
        do i = 1, nnew
            if (i /= gdof(ii)) then
                nout = nout + 1
                indices(nout) = i
            else
                ii = min(ii + 1, nbc)
            end if
        end do
    end if

    ! Remap the existing nonzeros without creating entries in constrained rows.
    allocate(rows(nnz), cols(nnz), vals(nnz))
    pos = 0
    do i = 1, n
        do j = x%row_indices(i), x%row_indices(i + 1) - 1
            pos = pos + 1
            rows(pos) = indices(i)
            cols(pos) = indices(x%column_indices(j))
            vals(pos) = x%values(j)
        end do
    end do
    rst = create_csr_matrix(nnew, nnew, rows, cols, vals)
end function

! ------------------------------------------------------------------------------
! REF: https://www.sciencedirect.com/topics/engineering/prescribed-displacement-boundary-condition
subroutine apply_displacement_constraint_dense(dof, val, k, f)
    !! Applies a displacement constraint to the specified degree of freedom.
    integer(int32), intent(in) :: dof
        !! The global degree-of-freedom to which the constraint should be
        !! applied.
    real(real64), intent(in) :: val
        !! The value of the displacement constraint.
    real(real64), intent(inout), dimension(:,:) :: k
        !! The stiffness matrix to which the constraint should be applied.
    real(real64), intent(inout), dimension(:) :: f
        !! The external force vector to which the constraint should be applied.

    ! Wipe out the rows in the matrix and place a value of 1 on the diagonal
    k(dof,:) = 0.0d0
    k(dof,dof) = 1.0d0

    ! Update the external force vector
    f(dof) = val
end subroutine

! ------------------------------------------------------------------------------
subroutine apply_displacement_constraint_csr(dof, val, k, f)
    !! Applies a displacement constraint to a CSR-format sparse matrix.
    integer(int32), intent(in) :: dof
        !! The global degree-of-freedom to which the constraint should be
        !! applied.
    real(real64), intent(in) :: val
        !! The value of the displacement constraint.
    type(csr_matrix), intent(inout) :: k
        !! The stiffness matrix to which the constraint should be applied.
    real(real64), intent(inout), dimension(:) :: f
        !! The external force vector to which the constraint should be applied.

    ! Local Variables
    integer(int32) :: i, j, m, n, nnz, pos
    integer(int32), allocatable :: rows(:), cols(:)
    real(real64), allocatable :: vals(:)

    ! Initialization
    m = size(k, 1)
    n = size(k, 2)
    nnz = size(k%values)

    ! Input Checking
    if (m /= n) error stop DYN_MATRIX_SIZE_ERROR
    if (dof < 1 .or. dof > m) error stop DYN_INDEX_OUT_OF_RANGE
    if (size(f) /= m) error stop DYN_ARRAY_SIZE_ERROR

    ! Rebuild the sparse matrix, omitting the constrained row and inserting its
    ! unit diagonal entry.
    allocate(rows(nnz + 1), cols(nnz + 1), vals(nnz + 1))
    pos = 0
    do i = 1, m
        if (i == dof) cycle
        do j = k%row_indices(i), k%row_indices(i + 1) - 1
            pos = pos + 1
            rows(pos) = i
            cols(pos) = k%column_indices(j)
            vals(pos) = k%values(j)
        end do
    end do
    pos = pos + 1
    rows(pos) = dof
    cols(pos) = dof
    vals(pos) = 1.0d0
    k = create_csr_matrix(m, n, rows(1:pos), cols(1:pos), vals(1:pos))
    f(dof) = val
end subroutine

! ******************************************************************************
! ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function e_stiffness_matrix(this, rule) result(rst)
    !! Computes the stiffness matrix for the element.
    class(element), intent(in) :: this
        !! The element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: r

    ! Initialization
    if (present(rule)) then
        r = rule
    else
        r = DYN_TWO_POINT_INTEGRATION_RULE
    end if

    ! Process
    rst = integrate(element_stiffness_integrand, this, r)
end function

! ----------
pure function element_stiffness_integrand(elem, s) result(rst)
    !! The integrand function for computing the stiffness matrix of an element.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate vector at which to evaluate the integrand.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integrand.

    ! Local Variables
    real(real64) :: jdet
    real(real64), allocatable, dimension(:,:) :: b, bt, d, x, jac

    ! Process
    b = elem%strain_displacement_matrix(s)
    bt = transpose(b)
    d = elem%constitutive_matrix()
    jac = elem%jacobian(s)
    jdet = det(jac)
    x = matmul(d, b)
    rst = jdet * matmul(bt, x)
end function

! ------------------------------------------------------------------------------
pure function e_mass_matrix(this, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(element), intent(in) :: this
        !! The element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    integer(int32) :: r

    ! Initialization
    if (present(rule)) then
        r = rule
    else
        r = DYN_TWO_POINT_INTEGRATION_RULE
    end if

    ! Process
    rst = integrate(element_mass_integrand, this, r)
end function

! ----------
pure function element_mass_integrand(elem, s) result(rst)
    !! The integrand function for computing the mass matrix of an element.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate vector at which to evaluate the integrand.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integrand.

    ! Local Variables
    real(real64) :: jdet
    real(real64), allocatable, dimension(:,:) :: N, Nt, jac

    ! Process
    N = elem%shape_function_matrix(s)
    Nt = transpose(N)
    jac = elem%jacobian(s)
    jdet = det(jac)
    rst = elem%material%density * jdet * matmul(Nt, N)
end function

! ------------------------------------------------------------------------------
pure function e_ext_force_vector(this, q, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(element), intent(in) :: this
        !! The element object.
    real(real64), intent(in), dimension(:) :: q
        !! The surface traction forces vector or body force vector.  
        !! For instance, a 2D problem this vector would look like [qx, qy]**T.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    integer(int32) :: r

    ! Initialization
    if (present(rule)) then
    r = rule
    else
    r = DYN_TWO_POINT_INTEGRATION_RULE
    end if

    ! Process
    rst = matmul( &
        integrate(element_ext_force_integrand, this, r), &
        q &
    )
end function

! ----------
pure function element_ext_force_integrand(elem, s) result(rst)
    !! The integrand function for computing the external force vector of an 
    !! element.
    class(element), intent(in) :: elem
        !! The element object.
    real(real64), intent(in), dimension(:) :: s
        !! The natural coordinate vector at which to evaluate the integrand.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The integrand.

    ! Local Variables
    real(real64) :: jdet
    real(real64), allocatable, dimension(:,:) :: Nt, jac

    ! Process
    Nt = transpose(elem%shape_function_matrix(s))
    jac = elem%jacobian(s)
    jdet = det(jac)
    rst = jdet * Nt
end function

! ******************************************************************************
! LINE_ELEMENT MEMBERS
! ------------------------------------------------------------------------------
pure function le_length(this) result(rst)
    !! Computes the length of the line_element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    real(real64) :: rst
        !! The length of the line element.

    ! Local Variables
    real(real64) :: dx, dy, dz
    integer(int32) :: i1, i2
    type(node) :: n1, n2

    ! Process
    call this%get_terminal_nodes(i1, i2)
    n1 = this%get_node(i1)
    n2 = this%get_node(i2)
    dx = n2%x - n1%x
    dy = n2%y - n1%y
    dz = n2%z - n1%z
    rst = sqrt(dx**2 + dy**2 + dz**2)
end function

! ------------------------------------------------------------------------------
pure function le_stiffness_matrix(this, rule) result(rst)
    !! Computes the stiffness matrix for the element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Compute the stiffness matrix and apply the rotation transformation
    rst = e_stiffness_matrix(this, rule)
    rst = matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
pure function le_mass_matrix(this, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - MECH_ONE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_TWO_POINT_INTEGRATION_RULE
        !!
        !! - MECH_THREE_POINT_INTEGRATION_RULE
        !!
        !! - MECH_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is MECH_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting matrix.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: T, Tt

    ! Compute the rotation matrix
    T = this%rotation_matrix()
    Tt = transpose(T)

    ! Compute the mass matrix and apply the rotation transformation
    rst = e_mass_matrix(this, rule)
    rst = this%area * matmul(Tt, matmul(rst, T))
end function

! ------------------------------------------------------------------------------
pure function le_ext_force_vector(this, q, rule) result(rst)
    !! Computes the mass matrix for the element.
    class(line_element), intent(in) :: this
        !! The line_element object.
    real(real64), intent(in), dimension(:) :: q
        !! The surface traction forces vector or body force vector.  
        !! For instance, a 2D problem this vector would look like [qx, qy]**T.
    integer(int32), intent(in), optional :: rule
        !! The integration rule.  The rule must be one of the following:
        !!
        !! - DYN_ONE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_TWO_POINT_INTEGRATION_RULE
        !!
        !! - DYN_THREE_POINT_INTEGRATION_RULE
        !!
        !! - DYN_FOUR_POINT_INTEGRATION_RULE
        !!
        !! The default integration rule is DYN_TWO_POINT_INTEGRATION_RULE.
    real(real64), allocatable, dimension(:) :: rst
        !! The resulting vector.

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: T

    ! Compute the rotation matrix
    T = this%rotation_matrix()

    ! Compute the force vector
    rst = e_ext_force_vector(this, q, rule)
    rst = matmul(T, rst)
end function

! ******************************************************************************
! NODE MEMBERS
! ------------------------------------------------------------------------------
pure function nd_init_1(index, dof, x, y, z) result(rst)
    !! Constructs a new [[node]].
    integer(int32), intent(in) :: index
        !! The global index of the node.
    integer(int32), intent(in) :: dof
        !! The number of degrees of freedom of the node.
    real(real64), intent(in) :: x
        !! The x-coordinate.
    real(real64), intent(in) :: y
        !! The y-coordinate.
    real(real64), intent(in) :: z
        !! The z-coordinate.
    type(node) :: rst
        !! The new [[node]].
    rst%index = index
    rst%dof = dof
    rst%x = x
    rst%y = y
    rst%z = z
end function

! ------------------------------------------------------------------------------
pure function nd_init_2(index, dof, pt) result(rst)
    !! Constructs a new [[node]].
    integer(int32), intent(in) :: index
        !! The global index of the node.
    integer(int32), intent(in) :: dof
        !! The number of degrees of freedom of the node.
    class(point), intent(in) :: pt
        !! The location of the node.
    type(node) :: rst
        !! The new [[node]].
    rst = nd_init_1(index, dof, pt%x, pt%y, pt%z)
end function

! ******************************************************************************
! MATERIAL MEMBERS
! ------------------------------------------------------------------------------
pure function mat_init(modulus, pratio, density) result(rst)
    !! Constructs a new [[material]].
    real(real64), intent(in) :: modulus
        !! The modulus of elasticity.
    real(real64), intent(in) :: pratio
        !! The Poisson's ratio.
    real(real64), intent(in) :: density
        !! The density.
    type(material) :: rst
        !! The new [[material]].
    rst%modulus = modulus
    rst%poissons_ratio = pratio
    rst%density = density
end function



! ******************************************************************************
! PRIVATE ROUTINES
! ------------------------------------------------------------------------------
pure function det_1(x) result(rst)
    ! Determinant of a 1-by-1 matrix.
    real(real64), intent(in), dimension(:,:) :: x
    real(real64) :: rst
    rst = x(1,1)
end function

! ------------------------------------------------------------------------------
pure function det_2(x) result(rst)
    ! Determinant of a 2-by-2 matrix.
    real(real64), intent(in), dimension(:,:) :: x
    real(real64) :: rst

    rst = x(1,1) * x(2,2) - x(1,2) * x(2,1)
end function

! ------------------------------------------------------------------------------
pure function det_3(x) result(rst)
    ! Determinant of a 3-by-3 matrix.
    real(real64), intent(in), dimension(:,:) :: x
    real(real64) :: rst
    rst = x(1,1) * (x(2,2) * x(3,3) - x(2,3) * x(3,2)) - &
        x(1,2) * (x(2,1) * x(3,3) - x(2,3) * x(3,1)) + &
        x(1,3) * (x(2,1) * x(3,2) - x(2,2) * x(3,1))
end function

! ------------------------------------------------------------------------------
pure function det(x) result(rst)
    !! Computes the determinant of a matrix.
    real(real64), intent(in), dimension(:,:) :: x
        !! The matrix on which to operate.
    real(real64) :: rst
        !! The determinant.

    select case (size(x, 1))
    case (1)
        rst = det_1(x)
    case (2)
        rst = det_2(x)
    case (3)
        rst = det_3(x)
    case default
        rst = 0.0d0
    end select
end function

! ******************************************************************************
! SOLVERS
! ------------------------------------------------------------------------------
pure function solve_static_system_dense(K, F) result(rst)
    !! Solves the static system \(K u = f\).
    real(real64), intent(in), dimension(:,:) :: K
        !! The N-by-N stiffness matrix.
    real(real64), intent(in), dimension(:) :: F
        !! The N-element external forcing vector.
    real(real64), allocatable, dimension(:) :: rst
        !! The N-element solution vector.

    ! Local Variables
    integer(int32) :: n
    integer(int32), allocatable, dimension(:) :: pvt
    real(real64), allocatable, dimension(:,:) :: lu

    ! Input Check
    n = size(K, 1)
    if (size(K, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
    if (size(F) /= n) error stop DYN_ARRAY_SIZE_ERROR

    ! Factor the system
    call lu_factor(K, ipvt = pvt, lu = lu)

    ! Solve the system
    rst = solve_lu(lu, pvt, F)
end function

! ------------------------------------------------------------------------------
pure function solve_static_system_csr(K, F) result(rst)
    !! Solves the static system \(K u = f\).
    type(csr_matrix), intent(in) :: K
        !! The N-by-N stiffness matrix.
    real(real64), intent(in), dimension(:) :: F
        !! The N-element external forcing vector.
    real(real64), allocatable, dimension(:) :: rst
        !! The N-element solution vector.

    ! Local Variables
    integer(int32) :: n

    ! Input Check
    n = size(K, 1)
    if (size(K, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
    if (size(F) /= n) error stop DYN_ARRAY_SIZE_ERROR

    ! Solve the system
    rst= sparse_direct_solve(K, F)
end function

! ------------------------------------------------------------------------------
end module