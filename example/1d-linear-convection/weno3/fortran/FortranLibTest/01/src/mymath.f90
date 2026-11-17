! 简单的数学库：提供加法和乘法函数
module mymath_module
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: add, multiply

contains

    ! 实数加法
    function add(a, b) result(res)
        real(real64), intent(in) :: a, b
        real(real64) :: res
        res = a + b
    end function add

    ! 实数乘法
    function multiply(a, b) result(res)
        real(real64), intent(in) :: a, b
        real(real64) :: res
        res = a * b
    end function multiply

end module mymath_module