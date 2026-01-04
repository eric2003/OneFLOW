! 测试程序：调用mymath库的函数，验证功能
program test_mymath
    use iso_fortran_env, only: real64  ! 显式引入，避免模块依赖缺失导致的real64未定义
    use mymath_module
    implicit none
    real(real64) :: a, b, sum_res, mul_res

    ! 测试数据
    a = 2.5_real64
    b = 4.0_real64

    ! 调用库函数
    sum_res = add(a, b)
    mul_res = multiply(a, b)

    ! 输出结果
    print *, "=== Testing mymath library ==="
    print *, "a = ", a
    print *, "b = ", b
    print *, "a + b = ", sum_res
    print *, "a * b = ", mul_res
    print *, "=== Test completed ==="

    ! 简单验证（确保结果正确）
    if (abs(sum_res - 6.5_real64) < 1e-10_real64 .and. &
        abs(mul_res - 10.0_real64) < 1e-10_real64) then
        print *, "✓ Test passed!"
    else
        print *, "✗ Test failed!"
        stop 1
    end if

end program test_mymath