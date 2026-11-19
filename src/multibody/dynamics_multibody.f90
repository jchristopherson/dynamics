module dynamics_multibody
    use iso_fortran_env
    use dynamics_flexible_element
    use dynamics_linkage
    use dynamics_parallel_linkage
    use dynamics_joints
    implicit none
    private

    type :: multibody_model
        class(link), allocatable, dimension(:) :: m_rigidLinks
            !! A collection of rigid links that make up the model.
        class(flexible_element), allocatable, dimension(:) :: m_elements
            !! A collection of the flexible elements that make up the model.
        type(joint), allocatable, dimension(:) :: m_joints
            !! A collection of joints.
        real(real64), allocatable, dimension(:,:,:) :: m_parentFrames
            !! The parent coordinate frames.
        real(real64), allocatable, dimension(:,:,:) :: m_childFrames
            !! The child coordinate frames.
    contains
        ! TO DO:
        ! - assemble global mass matrix
        ! - assemble global internal force vector
        ! - assemble global external force vector
        ! - allow addition of rigid links, elements, joints, etc.
    end type
end module