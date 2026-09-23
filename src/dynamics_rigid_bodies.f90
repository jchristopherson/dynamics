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
module dynamics_rigid_bodies
    use iso_fortran_env
    implicit none
    private
    public :: rigid_body
    public :: initialize_rigid_body

    type rigid_body
        !! Defines a rigid body with mass properties about its center of mass.
        !! The translational and rotational kinetic energies are represented by
        !! $$ T=\frac{1}{2}m\boldsymbol{v}_{cg}^{T}\boldsymbol{v}_{cg}
        !! +\frac{1}{2}\boldsymbol{\omega}^{T}I_{cg}\boldsymbol{\omega}. $$
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
    !! Initializes a rigid_body object. The inertia tensor is interpreted about
    !! the supplied center of gravity, not about the world origin.
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

    call initialize_rigid_body(rst, m, inertia, cg)
end function

! ------------------------------------------------------------------------------
pure subroutine initialize_rigid_body(bdy, m, inertia, cg)
    !! Initializes a rigid_body object.
    class(rigid_body), intent(inout) :: bdy
        !! The rigid_body object.
    real(real64), intent(in), optional :: m
        !! The mass of the body.  If no mass is specified, a value of 1 is used.
    real(real64), intent(in), optional :: inertia(3, 3)
        !! The 3-by-3 inertia tensor.  If not supplied, an identity matrix
        !! is used.
    real(real64), intent(in), optional :: cg(3)
        !! The x-y-z location of the CG relative to the body coordinate frame.
        !! If not supplied, the CG is set to (0, 0, 0).

    if (present(m)) then
        bdy%mass = m
    else
        bdy%mass = 1.0d0
    end if

    if (present(inertia)) then
        bdy%inertia = inertia
    else
        bdy%inertia = reshape( &
            [1.0d0, 0.0d0, 0.0d0, &
            0.0d0, 1.0d0, 0.0d0, &
            0.0d0, 0.0d0, 1.0d0], &
            [3, 3] &
        )
    end if

    if (present(cg)) then
        bdy%cg = cg
    else
        bdy%cg = 0.0d0
    end if
end subroutine

! ------------------------------------------------------------------------------
end module