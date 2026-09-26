module dynamics_c_quaternions
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
subroutine c_quaternion_from_array(x, q) &
    bind(C, name = "c_quaternion_from_array")
    real(c_double), intent(in) :: x(4)
    type(c_quaternion), intent(out) :: q
    q%w = x(1)
    q%x = x(2)
    q%y = x(3)
    q%z = x(4)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_from_matrix(x, ldx, q) &
    bind(C, name = "c_quaternion_from_matrix")
    integer(c_int), intent(in), value :: ldx
    real(c_double), intent(in) :: x(ldx,3)
    type(c_quaternion), intent(out) :: q

    q = quaternion(x(1:3,1:3))
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_from_angle_axis(angle, axis, q) &
    bind(C, name = "c_quaternion_from_angle_axis")
    real(c_double), intent(in), value :: angle
    real(c_double), intent(in) :: axis(3)
    type(c_quaternion), intent(out) :: q

    q = quaternion(angle, axis)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_normalize(q) bind(C, name = "c_quaternion_normalize")
    type(c_quaternion), intent(inout) :: q
    type(quaternion) :: qf
    qf = q
    call qf%normalize()
    q = qf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_add(x, y, q) bind(C, name = "c_quaternion_add")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf + yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_subtract(x, y, q) bind(C, name = "c_quaternion_subtract")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf - yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_multiply(x, y, q) bind(C, name = "c_quaternion_multiply")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf * yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_divide(x, y, q) bind(C, name = "c_quaternion_divide")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    q = xf / yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_scale(x, y, q) bind(C, name = "c_quaternion_scale")
    real(c_double), intent(in), value :: x
    type(c_quaternion), intent(in) :: y
    type(c_quaternion), intent(out) :: q
    type(quaternion) :: yf
    yf = y
    q = x * yf
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_conjugate(q, qc) bind(C, name = "c_quaternion_conjugate")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: qc
    type(quaternion) :: qf
    qf = q
    qc = conjg(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_rotate(q, r, rp) bind(C, name = "c_quaternion_rotate")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(in) :: r(3)
    real(c_double), intent(out) :: rp(3)
    type(quaternion) :: qf
    qf = q
    rp = aimag(qf * r * conjg(qf))
end subroutine

! ------------------------------------------------------------------------------
function c_quaternion_abs(q) result(rst) bind(C, name = "c_quaternion_abs")
    type(c_quaternion), intent(in) :: q
    real(c_double) :: rst
    type(quaternion) :: qf
    qf = q
    rst = abs(qf)
end function

! ------------------------------------------------------------------------------
subroutine c_quaternion_inverse(q, qinv) bind(C, name = "c_quaternion_inverse")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: qinv
    type(quaternion) :: qf
    qf = q
    qinv = inverse(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_to_matrix(q, r, ldr) bind(C, name = "c_quaternion_to_matrix") 
    type(c_quaternion), intent(in) :: q
    integer(c_int), intent(in), value :: ldr
    real(c_double), intent(out) :: r(ldr,3)
    type(quaternion) :: qf
    qf = q
    r(1:3,1:3) = qf%to_matrix()
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_to_angle_axis(q, angle, axis) &
    bind(C, name = "c_quaternion_to_angle_axis")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(out) :: angle
    real(c_double), intent(out) :: axis(3)
    type(quaternion) :: qf
    qf = q
    call qf%to_angle_axis(angle, axis)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_exp(q, rst) bind(C, name = "c_quaternion_exp")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: rst
    type(quaternion) :: qf
    qf = q
    rst = exp(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_log(q, rst)  bind(C, name = "c_quaternion_log")
    type(c_quaternion), intent(in) :: q
    type(c_quaternion), intent(out) :: rst
    type(quaternion) :: qf
    qf = q
    rst = log(qf)
end subroutine

! ------------------------------------------------------------------------------
subroutine c_quaternion_pow(q, exponent, rst) bind(C, name = "c_quaternion_pow")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(in), value :: exponent
    type(c_quaternion), intent(out) :: rst
    type(quaternion) :: qf
    qf = q
    rst = qf**exponent
end subroutine

! ------------------------------------------------------------------------------
function c_quaternion_dot_product(x, y) result(rst) &
    bind(C, name = "c_quaternion_dot_product")
    type(c_quaternion), intent(in) :: x
    type(c_quaternion), intent(in) :: y
    real(c_double) :: rst
    type(quaternion) :: xf, yf
    xf = x
    yf = y
    rst = dot_product(xf, yf)
end function

! ------------------------------------------------------------------------------
subroutine c_quaternion_to_roll_pitch_yaw(q, roll, pitch, yaw) &
    bind(C, name = "c_quaternion_to_roll_pitch_yaw")
    type(c_quaternion), intent(in) :: q
    real(c_double), intent(out) :: roll
    real(c_double), intent(out) :: pitch
    real(c_double), intent(out) :: yaw
    type(quaternion) :: qf
    qf = q
    call qf%to_roll_pitch_yaw(roll, pitch, yaw)
end subroutine

end module
