module dynamics_c_maps
    use iso_c_binding
    use iso_fortran_env
    use dynamics
    use diffeq
    use spectrum, only : window
    use dynamics_error_handling
    use nonlin
    use dynamics_c_types
    implicit none

    ! The active C coordinate callback; the Fortran interface has no user data.
    procedure(c_poincare_coordinates), pointer, private :: &
        active_coordinates => null()
    !$omp threadprivate(active_coordinates)

contains

subroutine c_poincare_map(n, x, y, z, pln, side, nbuffer, xbuff, ybuff, zbuff, &
    nactual) bind(C, name = "c_poincare_map")
    use dynamics_maps
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(in) :: z(n)
    type(c_plane), intent(in) :: pln
    integer(c_int), intent(in), value :: side
    integer(c_int), intent(in), value :: nbuffer
    real(c_double), intent(out) :: xbuff(nbuffer)
    real(c_double), intent(out) :: ybuff(nbuffer)
    real(c_double), intent(out) :: zbuff(nbuffer)
    integer(c_int), intent(out) :: nactual

    ! Local Variables
    integer(int32) :: i
    type(plane) :: p
    real(real64), allocatable, dimension(:,:) :: rst

    ! Process
    p = pln
    rst = poincare_map(x, y, z, p, side)
    nactual = min(nbuffer, size(rst, 1))
    do concurrent (i = 1:nactual)
        xbuff(i) = rst(i,1)
        ybuff(i) = rst(i,2)
        zbuff(i) = rst(i,3)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine c_poincare_map_ode(fcn, tspan, n, iv, sample_count, pln, side, &
    solver, coordinates, nbuffer, xbuff, ybuff, zbuff, nactual) &
    bind(C, name = "c_poincare_map_ode")
    type(c_funptr), intent(in), value :: fcn
    real(c_double), intent(in) :: tspan(2)
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: iv(n)
    integer(c_int), intent(in), value :: sample_count
    type(c_plane), intent(in) :: pln
    integer(c_int), intent(in), value :: side
    integer(c_int), intent(in), value :: solver
    type(c_funptr), intent(in), value :: coordinates
    integer(c_int), intent(in), value :: nbuffer
    real(c_double), intent(out) :: xbuff(nbuffer)
    real(c_double), intent(out) :: ybuff(nbuffer)
    real(c_double), intent(out) :: zbuff(nbuffer)
    integer(c_int), intent(out) :: nactual

    ! Local Variables
    integer(int32) :: i
    type(plane) :: p
    type(ode_container) :: sys
    type(c_ode_equations_container) :: arg
    procedure(c_ode_equations), pointer :: fptr
    type(runge_kutta_23), target :: rk23
    type(runge_kutta_45), target :: rk45
    type(runge_kutta_853), target :: rk853
    type(rosenbrock), target :: rbrk
    type(bdf), target :: bdiff
    type(adams), target :: pece
    type(kennedy_carpenter_4), target :: kc4
    type(kennedy_carpenter_5), target :: kc5
    type(tsitouras_54), target :: t54
    class(ode_integrator), pointer :: integrator_obj
    real(real64), allocatable, dimension(:,:) :: rst

    ! Initialization
    if (.not.c_associated(fcn)) error stop DYN_NULL_POINTER_ERROR
    call c_f_procpointer(fcn, fptr)
    arg%fcn => fptr
    sys%fcn => cpm_ode_fcn
    p = pln

    select case (solver)
    case (DYN_ADAMS)
        integrator_obj => pece
    case (DYN_BDF)
        integrator_obj => bdiff
    case (DYN_ROSENBROCK)
        integrator_obj => rbrk
    case (DYN_RUNGE_KUTTA_23)
        integrator_obj => rk23
    case (DYN_RUNGE_KUTTA_45)
        integrator_obj => rk45
    case (DYN_RUNGE_KUTTA_853)
        integrator_obj => rk853
    case (DYN_KENNEDY_CARPENTER_4)
        integrator_obj => kc4
    case (DYN_KENNEDY_CARPENTER_5)
        integrator_obj => kc5
    case (DYN_TSITOURAS_5)
        integrator_obj => t54
    case default
        integrator_obj => rk45
    end select

    ! Process
    if (c_associated(coordinates)) then
        call c_f_procpointer(coordinates, active_coordinates)
        rst = poincare_map(sys, tspan, iv, sample_count, pln = p, &
            side = side, solver = integrator_obj, &
            coordinates = cpm_coordinates, args = arg)
        active_coordinates => null()
    else
        rst = poincare_map(sys, tspan, iv, sample_count, pln = p, &
            side = side, solver = integrator_obj, args = arg)
    end if
    nactual = min(nbuffer, size(rst, 1))
    do concurrent (i = 1:nactual)
        xbuff(i) = rst(i,1)
        ybuff(i) = rst(i,2)
        zbuff(i) = rst(i,3)
    end do
end subroutine

! --------------------
subroutine cpm_ode_fcn(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args
    select type (args)
    class is (c_ode_equations_container)
        call args%fcn(size(y), x, y, dydx)
    end select
end subroutine

! --------------------
subroutine cpm_coordinates(t, state, coordinates_out)
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: state
    real(real64), intent(out), dimension(3) :: coordinates_out
    call active_coordinates(size(state), t, state, coordinates_out)
end subroutine

end module
