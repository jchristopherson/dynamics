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
module dynamics_flexible_joints
    !! Kinematic joints between node-table positions of a flexible mesh.
    !!
    !! A joint couples two nodes through reference frames attached to each
    !! node. The current joint frame on node `i` is `exp([theta_i]) * F_i`,
    !! where `theta_i` is the nodal rotation vector and `F_i` is the reference
    !! frame. This matches the rotation convention of the beam elements.
    !! Constraints are written in the current configuration and returned with
    !! their derivatives with respect to the four attachment vectors
    !! `[u1, theta1, u2, theta2]`.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dynamics_error_handling, only : DYN_ARRAY_SIZE_ERROR, DYN_INVALID_INPUT_ERROR
    use dynamics_joints, only : FIXED_JOINT, REVOLUTE_JOINT, PRISMATIC_JOINT, &
        CYLINDRICAL_JOINT, UNIVERSAL_JOINT, SPHERICAL_JOINT
    implicit none
    private

    public :: flexible_joint
    public :: joint_constraint_count
    public :: flexible_joint_constraints

    real(real64), parameter, dimension(3,3) :: identity_frame = reshape( &
        [1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3, 3])

    type :: flexible_joint
        !! Couples two nodes of a flexible mesh. Each frame is the reference
        !! orientation of the joint on its node, with columns as the joint
        !! x, y, and z axes.
        integer :: joint_type = REVOLUTE_JOINT
            !! One of the joint type constants from [[dynamics_joints]].
        integer :: node1 = 0
            !! Node-table position of the first attachment.
        integer :: node2 = 0
            !! Node-table position of the second attachment.
        real(real64), dimension(3,3) :: frame1 = identity_frame
            !! Reference joint frame on `node1`.
        real(real64), dimension(3,3) :: frame2 = identity_frame
            !! Reference joint frame on `node2`.
    end type

    interface flexible_joint
        module procedure :: joint_init
    end interface

contains
! ------------------------------------------------------------------------------
    function joint_init(joint_type, node1, node2, frame1, frame2) result(rst)
        !! Construct a validated joint. Omitted frames are aligned with the
        !! global axes. Invalid types, identical nodes, and frames that are
        !! not right-handed orthonormal cause an error stop.
        integer, intent(in) :: joint_type
            !! Joint type constant from [[dynamics_joints]].
        integer, intent(in) :: node1
            !! Node-table position of the first attachment.
        integer, intent(in) :: node2
            !! Node-table position of the second attachment.
        real(real64), intent(in), optional, dimension(3,3) :: frame1
            !! Reference joint frame on `node1`.
        real(real64), intent(in), optional, dimension(3,3) :: frame2
            !! Reference joint frame on `node2`.
        type(flexible_joint) :: rst
            !! The new joint.

        if (joint_constraint_count(joint_type) == 0) error stop DYN_INVALID_INPUT_ERROR
        if (node1 < 1 .or. node2 < 1 .or. node1 == node2) error stop DYN_INVALID_INPUT_ERROR
        rst%joint_type = joint_type
        rst%node1 = node1
        rst%node2 = node2
        if (present(frame1)) rst%frame1 = frame1
        if (present(frame2)) rst%frame2 = frame2
        if (.not.is_rotation(rst%frame1) .or. .not.is_rotation(rst%frame2)) then
            error stop DYN_INVALID_INPUT_ERROR
        end if
    end function

! ------------------------------------------------------------------------------
    pure function joint_constraint_count(joint_type) result(rst)
        !! Return the number of scalar constraints removed by a joint:
        !! six minus its relative degrees of freedom. Unknown types return zero.
        integer, intent(in) :: joint_type
            !! Joint type constant from [[dynamics_joints]].
        integer :: rst
            !! Number of constraint rows.

        select case (joint_type)
        case (FIXED_JOINT)
            rst = 6
        case (REVOLUTE_JOINT, PRISMATIC_JOINT)
            rst = 5
        case (CYLINDRICAL_JOINT, UNIVERSAL_JOINT)
            rst = 4
        case (SPHERICAL_JOINT)
            rst = 3
        case default
            rst = 0
        end select
    end function

! ------------------------------------------------------------------------------
    subroutine flexible_joint_constraints(joint, dx, theta1, theta2, g, jac)
        !! Evaluate the constraint values and their derivatives.
        !!
        !! The point difference `dx` is `(x1 + u1) - (x2 + u2)`, where `x_i`
        !! are reference positions and `u_i` are nodal displacements. The
        !! rotation vectors `theta_i` are the nodal rotation DOFs. Columns of
        !! `jac` are ordered `[u1, theta1, u2, theta2]`.
        !!
        !! Constraint rows are:
        !!
        !! - fixed: point coincidence and orientation equality
        !! - revolute: point coincidence and parallel z axes
        !! - prismatic: zero offset perpendicular to the joint z axis and
        !!   orientation equality
        !! - cylindrical: zero perpendicular offset and parallel z axes
        !! - universal: point coincidence and perpendicular z axes
        !! - spherical: point coincidence
        !!
        !! The orientation rows are singular at a relative rotation of pi.
        type(flexible_joint), intent(in) :: joint
            !! Joint definition with validated reference frames.
        real(real64), intent(in), dimension(3) :: dx
            !! Current point difference between the two attachments.
        real(real64), intent(in), dimension(3) :: theta1
            !! Rotation vector of the node on `joint%node1`.
        real(real64), intent(in), dimension(3) :: theta2
            !! Rotation vector of the node on `joint%node2`.
        real(real64), intent(out), dimension(:) :: g
            !! Constraint values, of length [[joint_constraint_count]].
        real(real64), intent(out), dimension(:,:) :: jac
            !! Jacobian with one row per constraint and 12 columns.

        real(real64), dimension(6,3) :: gu1
        real(real64), dimension(6,3) :: gu2
        real(real64), dimension(6,3) :: gw1
        real(real64), dimension(6,3) :: gw2
        real(real64), dimension(3,3) :: axes1
        real(real64), dimension(3,3) :: axes2
        real(real64), dimension(3,3) :: l1
        real(real64), dimension(3,3) :: l2
        integer :: m

        m = joint_constraint_count(joint%joint_type)
        if (size(g) /= m .or. size(jac,1) /= m .or. size(jac,2) /= 12) error stop DYN_ARRAY_SIZE_ERROR

        g = 0.0d0
        gu1 = 0.0d0
        gu2 = 0.0d0
        gw1 = 0.0d0
        gw2 = 0.0d0
        axes1 = matmul(rotation_matrix(theta1), joint%frame1)
        axes2 = matmul(rotation_matrix(theta2), joint%frame2)
        l1 = left_jacobian(theta1)
        l2 = left_jacobian(theta2)

        select case (joint%joint_type)
        case (FIXED_JOINT)
            call point_rows(1)
            call orientation_rows(4)
        case (REVOLUTE_JOINT)
            call point_rows(1)
            call alignment_rows(4)
        case (PRISMATIC_JOINT)
            call offset_rows(1)
            call orientation_rows(3)
        case (CYLINDRICAL_JOINT)
            call offset_rows(1)
            call alignment_rows(3)
        case (UNIVERSAL_JOINT)
            call point_rows(1)
            call perpendicular_row(4)
        case (SPHERICAL_JOINT)
            call point_rows(1)
        end select

        ! Rotation-vector derivatives use omega = J_l(theta) * dtheta.
        jac(:,1:3) = gu1(1:m,:)
        jac(:,4:6) = matmul(gw1(1:m,:), l1)
        jac(:,7:9) = gu2(1:m,:)
        jac(:,10:12) = matmul(gw2(1:m,:), l2)
    contains
        subroutine point_rows(first)
            !! Coincidence of the two attachment points.
            integer, intent(in) :: first
                !! First row of the block.

            g(first:first+2) = dx
            gu1(first:first+2,:) = identity_frame
            gu2(first:first+2,:) = -identity_frame
        end subroutine

        subroutine offset_rows(first)
            !! Zero offset of the point difference along the first two axes.
            integer, intent(in) :: first
                !! First row of the block.

            real(real64), dimension(3) :: p
            integer :: k

            do k = 1, 2
                p = axes1(:,k)
                g(first+k-1) = dot_product(p, dx)
                gu1(first+k-1,:) = p
                gu2(first+k-1,:) = -p
                gw1(first+k-1,:) = cross3(p, dx)
            end do
        end subroutine

        subroutine alignment_rows(first)
            !! Parallel z axes, written as zero components of a3 x b3 along
            !! the first two axes of the first joint frame.
            integer, intent(in) :: first
                !! First row of the block.

            real(real64), dimension(3) :: c
            real(real64), dimension(3) :: p
            real(real64) :: s
            integer :: k

            c = cross3(axes1(:,3), axes2(:,3))
            s = dot_product(axes1(:,3), axes2(:,3))
            do k = 1, 2
                p = axes1(:,k)
                g(first+k-1) = dot_product(p, c)
                gw1(first+k-1,:) = cross3(p, c) + dot_product(p, axes1(:,3))*axes2(:,3) - s*p
                gw2(first+k-1,:) = s*p - dot_product(p, axes2(:,3))*axes1(:,3)
            end do
        end subroutine

        subroutine perpendicular_row(row)
            !! Perpendicular z axes, written as a3 . b3 = 0.
            integer, intent(in) :: row
                !! Row of the constraint.

            real(real64), dimension(3) :: c

            c = cross3(axes1(:,3), axes2(:,3))
            g(row) = dot_product(axes1(:,3), axes2(:,3))
            gw1(row,:) = c
            gw2(row,:) = -c
        end subroutine

        subroutine orientation_rows(first)
            !! Equal joint frames, written as E = 1/2 sum(a_k x b_k).
            integer, intent(in) :: first
                !! First row of the block.

            real(real64), dimension(3) :: a
            real(real64), dimension(3) :: b
            real(real64) :: s
            integer :: k

            do k = 1, 3
                a = axes1(:,k)
                b = axes2(:,k)
                s = dot_product(a, b)
                g(first:first+2) = g(first:first+2) + 0.5d0*cross3(a, b)
                gw1(first:first+2,:) = gw1(first:first+2,:) + &
                    0.5d0*(spread(a, 2, 3)*spread(b, 1, 3) - s*identity_frame)
                gw2(first:first+2,:) = gw2(first:first+2,:) + &
                    0.5d0*(s*identity_frame - spread(b, 2, 3)*spread(a, 1, 3))
            end do
        end subroutine
    end subroutine

! ------------------------------------------------------------------------------
    pure function is_rotation(r) result(rst)
        !! True if `r` is finite, orthonormal, and right-handed.
        real(real64), intent(in), dimension(3,3) :: r
            !! Candidate frame.
        logical :: rst
            !! Validation result.

        rst = all(ieee_is_finite(r))
        if (.not.rst) return
        rst = maxval(abs(matmul(transpose(r), r) - identity_frame)) <= 1.0d-10
        rst = rst .and. abs(dot_product(cross3(r(:,1), r(:,2)), r(:,3)) - 1.0d0) <= 1.0d-10
    end function

! ------------------------------------------------------------------------------
    pure function cross3(a, b) result(c)
        !! Cross product of two three-vectors.
        real(real64), intent(in), dimension(3) :: a
            !! Left operand.
        real(real64), intent(in), dimension(3) :: b
            !! Right operand.
        real(real64), dimension(3) :: c
            !! Cross product.

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function

! ------------------------------------------------------------------------------
    pure function skew(v) result(k)
        !! Skew-symmetric matrix such that `skew(v) * x = v x x`.
        real(real64), intent(in), dimension(3) :: v
            !! Vector to map.
        real(real64), dimension(3,3) :: k
            !! Skew-symmetric matrix.

        k = reshape([0.0d0, v(3), -v(2), -v(3), 0.0d0, v(1), v(2), -v(1), 0.0d0], [3, 3])
    end function

! ------------------------------------------------------------------------------
    pure function rotation_matrix(theta) result(r)
        !! Rodrigues formula for the rotation vector `theta`.
        real(real64), intent(in), dimension(3) :: theta
            !! Rotation vector in radians.
        real(real64), dimension(3,3) :: r
            !! Rotation matrix.

        real(real64), dimension(3,3) :: k
        real(real64) :: phi
        real(real64) :: a
        real(real64) :: b

        phi = norm2(theta)
        k = skew(theta)
        if (phi < 1.0d-8) then
            a = 1.0d0
            b = 0.5d0
        else
            a = sin(phi)/phi
            b = 2.0d0*sin(0.5d0*phi)**2/phi**2
        end if
        r = identity_frame + a*k + b*matmul(k, k)
    end function

! ------------------------------------------------------------------------------
    pure function left_jacobian(theta) result(l)
        !! Left Jacobian of SO(3), so that a perturbation dtheta gives a
        !! spatial angular increment `J_l(theta) * dtheta`.
        real(real64), intent(in), dimension(3) :: theta
            !! Rotation vector in radians.
        real(real64), dimension(3,3) :: l
            !! Left Jacobian.

        real(real64), dimension(3,3) :: k
        real(real64) :: phi
        real(real64) :: phi2
        real(real64) :: a
        real(real64) :: b

        phi = norm2(theta)
        phi2 = phi*phi
        k = skew(theta)
        if (phi < 1.0d-2) then
            ! Series avoid cancellation in (phi - sin(phi)) / phi^3.
            a = 0.5d0 - phi2/24.0d0 + phi2*phi2/720.0d0
            b = 1.0d0/6.0d0 - phi2/120.0d0 + phi2*phi2/5040.0d0
        else
            a = (1.0d0 - cos(phi))/phi2
            b = (phi - sin(phi))/(phi2*phi)
        end if
        l = identity_frame + a*k + b*matmul(k, k)
    end function
end module
