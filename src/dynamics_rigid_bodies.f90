module dynamics_rigid_bodies
    use iso_fortran_env
    implicit none
    private
    public :: rigid_body

    type rigid_body
        !! Defines a rigid body.
        real(real64), public :: mass
            !! The mass of the body.
        real(real64), public :: cg(3)
            !! The x-y-z location of the CG relative to the body coordinate 
            !! frame.
        real(real64), public :: inertia(3, 3)
            !! The 3-by-3 inertia tensor as measured about the CG of the body.
    end type

    interface rigid_body
        module procedure :: rb_init
    end interface

contains
! ------------------------------------------------------------------------------
pure function rb_init(m, inertia, cg) result(rst)
    !! Initializes a rigid_body object.
    real(real64), intent(in), optional :: m
        !! The mass of the body.  If no mass is specified, a value of 1 is used.
    real(real64), intent(in), optional :: inertia(3, 3)
        !! The 3-by-3 inertia tensor.  If not supplied, an identity matrix
        !! is used.
    real(real64), intent(in), optional :: cg(3)
        !! The x-y-z location of the CG relative to the body coordinate frame.
        !! If not supplied, the CG is set to (0, 0, 0).
    type(rigid_body) :: rst
        !! The rigid_body object.

    if (present(m)) then
        rst%mass = m
    else
        rst%mass = 1.0d0
    end if

    if (present(inertia)) then
        rst%inertia = inertia
    else
        rst%inertia = reshape( &
            [1.0d0, 0.0d0, 0.0d0, &
            0.0d0, 1.0d0, 0.0d0, &
            0.0d0, 0.0d0, 1.0d0], &
            [3, 3] &
        )
    end if

    if (present(cg)) then
        rst%cg = cg
    else
        rst%cg = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
end module