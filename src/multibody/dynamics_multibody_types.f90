module dynamics_multibody_types
    use iso_fortran_env
    implicit none
    private
    public :: degree_of_freedom
    public :: multibody_node

    type :: degree_of_freedom
        !! Defines a single degree of freedom.
        integer :: equation_id
            !! The equation ID.
        integer :: node_id
            !! The associated multibody_node ID.
        integer :: component
            !! The component of the displacement vector.  For example:
            !! 1 = UX
            !! 2 = UY
            !! 3 = UZ
            !! 4 = RX
            !! 5 = RY
            !! 6 = RZ
    end type

    type :: multibody_node
        !! Defines a multibody node.
        integer :: id
            !! The node global ID number.
        real(real64), dimension(3) :: coordinates
            !! The [x, y, z] coordinates of the node.
    end type
end module