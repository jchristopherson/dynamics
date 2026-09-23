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
module example_models
    use iso_fortran_env
    use dynamics
    implicit none
contains
    ! 2nd Order Test Problem
    ! x" + 2 * z * wn * x' + wn**2 * x = f(t)
    pure function example_2nd_order_forcing(t) result(rst)
        ! Arguments
        real(real64), intent(in) :: t
        real(real64) :: rst

        ! Process
        rst = chirp(t, 1.0d2, 5.0d0, 1.0d0, 1.0d2)
    end function

    pure subroutine example_2nd_order(t, x, dxdt)
        ! Arguments
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt

        ! Model Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: z = 1.0d-1
        real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1

        ! Local Variables
        real(real64) :: f

        ! Process
        f = example_2nd_order_forcing(t)
        dxdt(1) = x(2)
        dxdt(2) = f - (2.0d0 * z * wn * x(2) + wn**2 * x(1))
    end subroutine
end module