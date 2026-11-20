module dynamics_flexible_body_types
    use iso_fortran_env
    implicit none
    private
    public :: degree_of_freedom
    public :: fe_node
    public :: DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ

    integer, parameter :: DOF_UX = 1
    integer, parameter :: DOF_UY = 2
    integer, parameter :: DOF_UZ = 3
    integer, parameter :: DOF_RX = 4
    integer, parameter :: DOF_RY = 5
    integer, parameter :: DOF_RZ = 6

    type :: degree_of_freedom
        !! Describes a numbered DOF in the current mesh layout.
        integer :: equation_id
            !! One-based global equation number.
        integer :: node_id
            !! The associated `fe_node%id`, not its node-table position.
        integer :: component
            !! Physical component label, such as `DOF_UY`, not a local slot.
    end type

    type :: fe_node
        !! Defines a node for a finite element body.
        integer :: id
            !! The node global ID number.
        real(real64), dimension(3) :: coordinates
            !! The [x, y, z] coordinates of the node.
    end type
end module