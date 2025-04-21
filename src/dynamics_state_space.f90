module dynamics_state_space
    use iso_fortran_env
    use diffeq
    implicit none
    private
    public :: state_space

    type state_space
        !! Defines a state-space representation of a dynamic system.  This
        !! implementation takes the form:
        !!
        !! $$ \dot{x}(t) = A(t) x(t) + B(t) u(t) $$
        !! $$ y(t) = C(t) x(t) + D(t) u(t) $$
        !!
        !! Where:
        !! 
        !! - /( t /) denotes time.
        !!
        !! - /( x(t) /) is the state vector.
        !!
        !! - /( u(t) /) is the input vector.
        !!
        !! - /( y(t) /) is the output vector.
        real(real64), allocatable, dimension(:,:) :: A
            !! The N-by-N dynamics matrix, where N is the number of state
            !! variables.
        real(real64), allocatable, dimension(:,:) :: B
            !! The N-by-M input matrix, where M is the number of inputs.
        real(real64), allocatable, dimension(:,:) :: C
            !! The P-by-N output matrix, where P is the number of outputs.
        real(real64), allocatable, dimension(:,:) :: D
            !! The P-by-M feedthrough matrix.
    end type

contains
end module