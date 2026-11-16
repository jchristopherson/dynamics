module dynamics_flexible_types
    use iso_fortran_env
    use autodiff
    implicit none
    private
    public :: flexible_element
    public :: finite_element

    type, abstract :: flexible_element
        !! The base type for all flexible element types.
    contains
        procedure(get_integer_property), deferred, public, pass :: get_dof_count
            !! Gets the number of degrees of freedom for this element.
    end type

    abstract interface
        pure function get_integer_property(this) result(rst)
            !! Gets an integer property value from a [[flexible_element]] based
            !! type.
            import flexible_element
            class(flexible_element), intent(in) :: this
                !! The [[flexible_element]] type.
            integer :: rst
                !! The property value.
        end function
    end interface

    type, abstract, extends(flexible_element) :: finite_element
    ! TO DO: Evaluate shape functions - use AD
    ! TO DO: Evaluate strain energy - use AD
    end type

end module