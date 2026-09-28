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
module dynamics_c_structural_solvers
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use dynamics_error_handling
    use dynamics_c_types
    implicit none

contains

function get_structural_integrator(obj) result(rst)
    ! Resolves an opaque structural integrator handle.
    type(c_ptr), intent(in), value :: obj
    class(structural_integrator), pointer :: rst

    type(c_structural_integrator_container), pointer :: cont

    rst => null()
    if (.not.c_associated(obj)) return
    call c_f_pointer(obj, cont)
    if (allocated(cont%item)) rst => cont%item
end function

! ------------------------------------------------------------------------------
function c_create_dense_generalized_alpha_integrator(n, m, ldm, c, ldc, k, &
    ldk, rho_infinity) result(rst) &
    bind(C, name = "c_create_dense_generalized_alpha_integrator")
    integer(c_int), intent(in), value :: n, ldm, ldc, ldk
    real(c_double), intent(in) :: m(ldm,n)
    real(c_double), intent(in) :: c(ldc,n)
    real(c_double), intent(in) :: k(ldk,n)
    real(c_double), intent(in), value :: rho_infinity
    type(c_ptr) :: rst

    type(c_structural_integrator_container), pointer :: cont
    type(dense_generalized_alpha_integrator) :: integrator

    if (ldm < n .or. ldc < n .or. ldk < n) error stop DYN_INVALID_INPUT_ERROR

    call integrator%initialize(m(1:n,1:n), c(1:n,1:n), k(1:n,1:n), &
        rho_infinity)
    allocate(cont)
    allocate(cont%item, source = integrator)
    rst = c_loc(cont)
end function

! ------------------------------------------------------------------------------
subroutine c_free_structural_integrator(obj) &
    bind(C, name = "c_free_structural_integrator")
    type(c_ptr), intent(in), value :: obj

    type(c_structural_integrator_container), pointer :: cont

    if (.not.c_associated(obj)) return
    call c_f_pointer(obj, cont)
    if (allocated(cont%item)) deallocate(cont%item)
    deallocate(cont)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_structural_integrator_step(obj, n, force_current, force_next, &
    dt, displacement, velocity, acceleration) &
    bind(C, name = "c_structural_integrator_step")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: force_current(n)
    real(c_double), intent(in) :: force_next(n)
    real(c_double), intent(in), value :: dt
    real(c_double), intent(inout) :: displacement(n)
    real(c_double), intent(inout) :: velocity(n)
    real(c_double), intent(inout) :: acceleration(n)

    class(structural_integrator), pointer :: integrator

    integrator => get_structural_integrator(obj)
    if (.not.associated(integrator)) error stop DYN_NULL_POINTER_ERROR
    call integrator%step(force_current, force_next, dt, displacement, &
        velocity, acceleration)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_structural_integrator_solve(obj, n, npts, forces, ldf, dt, &
    displacement, velocity, acceleration) &
    bind(C, name = "c_structural_integrator_solve")
    type(c_ptr), intent(in), value :: obj
    integer(c_int), intent(in), value :: n, npts, ldf
    real(c_double), intent(in) :: forces(ldf,npts)
    real(c_double), intent(in), value :: dt
    real(c_double), intent(inout) :: displacement(n)
    real(c_double), intent(inout) :: velocity(n)
    real(c_double), intent(inout) :: acceleration(n)

    class(structural_integrator), pointer :: integrator

    if (ldf < n) error stop DYN_INVALID_INPUT_ERROR
    integrator => get_structural_integrator(obj)
    if (.not.associated(integrator)) error stop DYN_NULL_POINTER_ERROR
    call integrator%solve(forces(1:n,:), dt, displacement, velocity, &
        acceleration)
end subroutine

! ------------------------------------------------------------------------------
end module
