module mymath_module
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: add, multiply

contains

    function add(a, b) result(res)
        real(real64), intent(in) :: a, b
        real(real64) :: res
        res = a + b
    end function add

    function multiply(a, b) result(res)
        real(real64), intent(in) :: a, b
        real(real64) :: res
        res = a * b
    end function multiply

end module mymath_module