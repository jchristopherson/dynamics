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

end module
