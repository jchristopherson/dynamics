module dynamics_flexible_element
    !! Common contracts for flexible elements and their assembled DOFs.
    use iso_fortran_env
    use dynamics_flexible_body_types, only : fe_node
    implicit none
    private
    public :: flexible_element

    type, abstract :: flexible_element
        !! Defines the DOF mapping, mass, internal force and tangent stiffness
        !! interface of a flexible element. Concrete types define their DOF
        !! ordering and force sign convention.
    contains
        procedure(get_dof_count), deferred, public, pass :: dof_count
            !! Total number of element DOFs.
        procedure(get_dof_indices), deferred, public, pass :: dof_indices
            !! Global equation indices in element DOF order.
        procedure(get_mass_matrix), deferred, public, pass :: mass_matrix
            !! Element mass matrix in element DOF order.
        procedure(get_internal_force), deferred, public, pass :: internal_force_vector
            !! Internal force evaluated at the supplied DOFs.
        procedure(get_stiffness_matrix), deferred, public, pass :: stiffness_matrix
            !! Tangent stiffness evaluated at the supplied DOFs.
        procedure :: add_damping_force
            !! Add a velocity-dependent internal force. Default is no damping.
    end type

    abstract interface
        pure function get_dof_count(this) result(rst)
            !! Gets the number of degrees of freedom for the flexible element.
            import flexible_element
            class(flexible_element), intent(in) :: this
                !! The [[flexible_element]] type.
            integer :: rst
                !! The number of degrees of freedom.
        end function

        pure function get_dof_indices(this) result(rst)
            !! Return an allocated array of global equation indices in the
            !! element's DOF order.
            import flexible_element
            class(flexible_element), intent(in) :: this
                !! The [[flexible_element]] type.
            integer, allocatable, dimension(:) :: rst
                !! Array with `dof_count()` entries.
        end function

        subroutine get_mass_matrix(this, m, nodes)
            !! Compute the element mass matrix. The caller supplies an output
            !! array with one row and column for each element DOF.
            use iso_fortran_env, only : real64
            import flexible_element, fe_node
            class(flexible_element), intent(inout) :: this
                !! The [[flexible_element]] type.
            real(real64), intent(out), dimension(:,:) :: m
                !! `dof_count()`-by-`dof_count()` mass matrix.
            type(fe_node), intent(in), dimension(:) :: nodes
                !! Caller-owned reference node table.
        end subroutine

        subroutine get_internal_force(this, q, f, nodes)
            !! Compute the internal force vector at the supplied DOFs.
            !! Concrete elements specify the sign convention.
            use iso_fortran_env, only : real64
            import flexible_element, fe_node
            class(flexible_element), intent(inout), target :: this
                !! The [[flexible_element]] type.
            real(real64), intent(in), dimension(:) :: q
                !! `dof_count()` element displacement coordinates.
            real(real64), intent(out), dimension(:) :: f
                !! `dof_count()` forces conjugate to the element DOFs.
            type(fe_node), intent(in), dimension(:) :: nodes
                !! Caller-owned reference node table.
        end subroutine

        subroutine get_stiffness_matrix(this, q, k, nodes)
            !! Compute the tangent stiffness matrix at the supplied DOFs.
            use iso_fortran_env, only : real64
            import flexible_element, fe_node
            class(flexible_element), intent(inout), target :: this
                !! The [[flexible_element]] type.
            real(real64), intent(in), dimension(:) :: q
                !! `dof_count()` element displacement coordinates.
            real(real64), intent(out), dimension(:,:) :: k
                !! `dof_count()`-by-`dof_count()` tangent stiffness matrix.
            type(fe_node), intent(in), target, dimension(:) :: nodes
                !! Caller-owned reference node table, borrowed during evaluation.
        end subroutine
    end interface

contains

    subroutine add_damping_force(this, q, v, f, nodes)
        !! Add the element's damping contribution to an existing internal
        !! force vector. The base implementation leaves it unchanged.
        class(flexible_element), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: q
            !! Local displacements.
        real(real64), intent(in), dimension(:) :: v
            !! Local coordinate rates.
        real(real64), intent(inout), dimension(:) :: f
            !! Local internal force; subtract the assembled vector in momentum balance.
        type(fe_node), intent(in), dimension(:) :: nodes
    end subroutine

end module