module dynamics_finite_element
    !! Finite-element topology contracts for flexible elements.
    use iso_fortran_env
    use dynamics_flexible_element
    implicit none
    private
    public :: finite_element

    type, abstract, extends(flexible_element) :: finite_element
        !! Extends the flexible-element interface with node count, DOFs per
        !! node, and spatial dimensionality.
    contains
        procedure(get_node_count), deferred, public, pass :: node_count
            !! Number of nodes in the element.
        procedure(get_node_connectivity), deferred :: node_connectivity
            !! One-based mesh node positions in local DOF block order.
        procedure(get_dof_per_node), deferred, public, pass :: dof_per_node
            !! Number of DOFs at each node.
        procedure(get_element_dimensionality), deferred, public, pass :: dimensionality
            !! Number of spatial coordinate dimensions.
    end type

    abstract interface
        pure function get_node_connectivity(this) result(indices)
            !! Return the node positions used by automatic mesh numbering.
            import finite_element
            class(finite_element), intent(in) :: this
                !! Element whose local vectors and matrices are grouped by node.
            integer, allocatable, dimension(:) :: indices
                !! `node_count()` distinct one-based positions in the mesh
                !! node table (not node IDs), in local DOF block order.
        end function

        pure function get_node_count(this) result(rst)
            !! Gets the number of nodes for this element.
            import finite_element
            class(finite_element), intent(in) :: this
                !! The [[finite_element]] object.
            integer :: rst
                !! The node count.
        end function

        pure function get_dof_per_node(this) result(rst)
            !! Gets the number of degrees of freedom per node for this element.
            import finite_element
            class(finite_element), intent(in) :: this
                !! The [[finite_element]] object.
            integer :: rst
                !! The degree-of-freedom count per node.
        end function

        pure function get_element_dimensionality(this) result(rst)
            !! Return the number of spatial coordinate dimensions.
            import finite_element
            class(finite_element), intent(in) :: this
                !! The [[finite_element]] object.
            integer :: rst
                !! The element dimensionality.
        end function
    end interface
end module