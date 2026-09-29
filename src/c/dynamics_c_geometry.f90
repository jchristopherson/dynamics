module dynamics_c_geometry
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

! ------------------------------------------------------------------------------
subroutine c_plane_normal(pln, nrm) bind(C, name = "c_plane_normal")
    type(c_plane), intent(in) :: pln
    real(c_double), intent(out) :: nrm(3)
    type(plane) :: p
    p = pln
    nrm = plane_normal(p)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_from_3_points(pt1, pt2, pt3, pln) &
    bind(C, name = "c_plane_from_3_points")
    real(c_double), intent(in) :: pt1(3)
    real(c_double), intent(in) :: pt2(3)
    real(c_double), intent(in) :: pt3(3)
    type(c_plane), intent(out) :: pln
    pln = plane(pt1, pt2, pt3)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_from_point_and_normal(pt, nrm, pln) &
    bind(C, name = "c_plane_from_point_and_normal")
    real(c_double), intent(in) :: pt(3)
    real(c_double), intent(in) :: nrm(3)
    type(c_plane), intent(out) :: pln
    pln = plane(pt, nrm)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plane_from_points(n, pts, ldp, pln) &
    bind(C, name = "c_plane_from_points")
    integer(c_int), intent(in), value :: n, ldp
    real(c_double), intent(in) :: pts(ldp,3)
    type(c_plane), intent(out) :: pln
    pln = plane(pts(1:n,:))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_flip_plane_normal(pln) bind(C, name = "c_flip_plane_normal")
    type(c_plane), intent(inout) :: pln
    type(plane) :: p
    p = pln
    call p%flip_normal()
    pln = p
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_2_points(pt1, pt2, ln) bind(C, name = "c_line_from_2_points")
    real(c_double), intent(in) :: pt1(3)
    real(c_double), intent(in) :: pt2(3)
    type(c_line), intent(out) :: ln
    ln = line(pt1, pt2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_2_planes(p1, p2, ln) bind(C, name = "c_line_from_2_planes")
    type(c_plane), intent(in) :: p1
    type(c_plane), intent(in) :: p2
    type(c_line), intent(out) :: ln
    type(plane) :: pln1, pln2
    pln1 = p1
    pln2 = p2
    ln = line(pln1, pln2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_points(n, pts, ldp, ln) bind(C, name = "c_line_from_points")
    integer(c_int), intent(in), value :: n, ldp
    real(c_double), intent(in) :: pts(ldp, 3)
    type(c_line), intent(out) :: ln
    ln = line(pts(1:n,:))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_evaluate_line_position(ln, t, x) bind(C, name = "c_evaluate_line_position")
    type(c_line), intent(in) :: ln
    real(c_double), intent(in), value :: t
    real(c_double), intent(out) :: x(3)
    type(line) :: l
    l = ln
    x = l%evaluate(t)
end subroutine

! ------------------------------------------------------------------------------
function c_is_parallel_vectors(n, x, y, tol) result(rst) &
    bind(C, name = "c_is_parallel_vectors")
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: x(n)
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    rst = logical(is_parallel(x, y, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_parallel_lines(x, y, tol) result(rst) &
    bind(C, name = "c_is_parallel_lines")
    type(c_line), intent(in) :: x
    type(c_line), intent(in) :: y
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(line) :: xf, yf
    xf = x
    yf = y
    rst = logical(is_parallel(xf, yf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_parallel_planes(x, y, tol) result(rst) &
    bind(C, name = "c_is_parallel_planes")
    type(c_plane), intent(in) :: x
    type(c_plane), intent(in) :: y
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(plane) :: xf, yf
    xf = x
    yf = y
    rst = logical(is_parallel(xf, yf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_point_on_plane(pt, pln, tol) result(rst) &
    bind(C, name = "c_is_point_on_plane")
    real(c_double), intent(in) :: pt(3)
    type(c_plane), intent(in) :: pln
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(plane) :: pf
    pf = pln
    rst = logical(is_point_on_plane(pt, pf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_is_point_on_line(pt, ln, tol) result(rst) &
    bind(C, name = "c_is_point_on_line")
    real(c_double), intent(in) :: pt(3)
    type(c_line), intent(in) :: ln
    real(c_double), intent(in), value :: tol
    logical(c_bool) :: rst
    type(line) :: lf
    lf = ln
    rst = logical(is_point_on_line(pt, lf, tol), c_bool)
end function

! ------------------------------------------------------------------------------
function c_nearest_point_on_line(pt, ln) result(rst) &
    bind(C, name = "c_nearest_point_on_line")
    real(c_double), intent(in) :: pt(3)
    type(c_line), intent(in) :: ln
    real(c_double) :: rst
    type(line) :: lf
    lf = ln
    rst = nearest_point_on_line(pt, lf)
end function

! ------------------------------------------------------------------------------
function c_point_to_line_distance(pt, ln) result(rst) &
    bind(C, name = "c_point_to_line_distance")
    real(c_double), intent(in) :: pt(3)
    type(c_line), intent(in) :: ln
    real(c_double) :: rst
    type(line) :: lf
    lf = ln
    rst = point_to_line_distance(pt, lf)
end function

! ------------------------------------------------------------------------------
function c_point_to_plane_distance(pt, pln) result(rst) &
    bind(C, name = "c_point_to_plane_distance")
    real(c_double), intent(in) :: pt(3)
    type(c_plane), intent(in) :: pln
    real(c_double) :: rst
    type(plane) :: pf
    pf = pln
    rst = point_to_plane_distance(pt, pf)
end function

! ------------------------------------------------------------------------------
subroutine c_vector_plane_projection(x, pln, px) &
    bind(C, name = "c_vector_plane_projection")
    real(c_double), intent(in) :: x(3)
    type(c_plane), intent(in) :: pln
    real(c_double), intent(out) :: px(3)
    type(plane) :: pf
    pf = pln
    px = vector_plane_projection(x, pf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_point_plane_projection(pt, pln, ppt) &
    bind(C, name = "c_point_plane_projection")
    real(c_double), intent(in) :: pt(3)
    type(c_plane), intent(in) :: pln
    real(c_double), intent(out) :: ppt(3)
    type(plane) :: pf
    pf = pln
    ppt = point_plane_projection(pt, pf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_2_points(pt1, pt2, ln) &
    bind(C, name = "c_plucker_line_from_2_points")
    real(c_double), intent(in) :: pt1(3)
    real(c_double), intent(in) :: pt2(3)
    type(c_plucker_line), intent(out) :: ln
    ln = plucker_line(pt1, pt2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_line(src, ln) &
    bind(C, name = "c_plucker_line_from_line")
    type(c_line), intent(in) :: src
    type(c_plucker_line), intent(out) :: ln
    type(line) :: fsrc
    fsrc = src
    ln = plucker_line(fsrc)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_2_planes(p1, p2, ln) &
    bind(C, name = "c_plucker_line_from_2_planes")
    type(c_plane), intent(in) :: p1
    type(c_plane), intent(in) :: p2
    type(c_plucker_line), intent(out) :: ln
    type(plane) :: pln1, pln2
    pln1 = p1
    pln2 = p2
    ln = plucker_line(pln1, pln2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_from_array(x, ln) &
    bind(C, name = "c_plucker_line_from_array")
    real(c_double), intent(in) :: x(6)
    type(c_plucker_line), intent(out) :: ln
    ln = plucker_line(x, nrm = .true.)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_mtx_mult(n, x, ldx, ln, y) &
    bind(C, name = "c_plucker_line_mtx_mult")
    integer(c_int), intent(in), value :: n, ldx
    real(c_double), intent(in) :: x(ldx,6)
    type(c_plucker_line), intent(in) :: ln
    real(c_double), intent(out) :: y(n)
    type(plucker_line) :: l
    l = ln
    y = matmul(x(1:n,:), l)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_plucker_line_to_array(ln, x) &
    bind(C, name = "c_plucker_line_to_array")
    type(c_plucker_line), intent(in) :: ln
    real(c_double), intent(out) :: x(6)
    x(1:3) = ln%u
    x(4:6) = ln%m
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_common_normal(ln1, ln2, ln) &
    bind(C, name = "c_line_common_normal")
    type(c_line), intent(in) :: ln1
    type(c_line), intent(in) :: ln2
    type(c_line), intent(out) :: ln

    type(line) :: f1, f2
    f1 = ln1
    f2 = ln2
    ln = line_common_normal(f1, f2)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_do_lines_intersect(ln1, ln2, intersect, t1, t2, tol) &
    bind(C, name = "c_do_lines_intersect")
    type(c_line), intent(in) :: ln1
    type(c_line), intent(in) :: ln2
    logical(c_bool), intent(out) :: intersect
    real(c_double), intent(out) :: t1
    real(c_double), intent(out) :: t2
    real(c_double), intent(in), value :: tol

    logical :: check
    type(line) :: f1, f2
    f1 = ln1
    f2 = ln2
    call do_lines_intersect(f1, f2, check, t1, t2, tol)
    intersect = logical(check, c_bool)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_line_from_point_and_vector(pt, v, ln) &
    bind(C, name = "c_line_from_point_and_vector")
    real(c_double), intent(in) :: pt(3)
    real(c_double), intent(in) :: v(3)
    type(c_line), intent(out) :: ln
    ln = line_from_point_and_vector(pt, v)
end subroutine

end module
