program test_precise_parameters
    ! 新增：引入iso_c_binding模块，用于关联Windows C API
    use, intrinsic :: iso_c_binding
    use kinds, only: dp
    use complex_wave_precise_module, only: PreciseParams, complex_wave_precise, init_precise_params, &
                                           get_precise_params, calculate_beta
    implicit none
    
    ! 声明Windows API：SetConsoleOutputCP（设置控制台输出代码页）
    interface
        subroutine SetConsoleOutputCP(cp) bind(C, name='SetConsoleOutputCP')
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: cp  ! 代码页编号（UTF-8对应65001）
        end subroutine SetConsoleOutputCP
    end interface
    
    integer :: i, n
    real(dp) :: x
    real(dp), allocatable :: x_vals(:), u_vals(:)
    type(PreciseParams) :: p
    
    ! 核心修改：调用API设置控制台输出编码为UTF-8（代码页65001）
    ! 仅在Windows平台（VS2022）生效，不影响程序其他功能
    call SetConsoleOutputCP(65001_c_int)
    
    ! 初始化精确参数
    call init_precise_params()
    p = get_precise_params()
    
    ! 打印参数详细信息
    print *, "=== Parameter Details ==="
    print *, "β calculation:"
    print *, "  log(2) = ", log(2.0_dp)
    print *, "  δ = ", p%delta
    print *, "  δ² = ", p%delta**2
    print *, "  36δ² = ", 36.0_dp * p%delta**2
    print *, "  β = log(2) / (36δ²) = ", p%beta
    print *, ""
    
    ! 验证β的计算
    print *, "=== Verification of β ==="
    print *, "For δ = 0.005:"
    print *, "  Manual calculation:"
    print *, "    δ² = 0.005² = 2.5e-5"
    print *, "    36δ² = 36 * 2.5e-5 = 9.0e-4"
    print *, "    log(2) ≈ 0.693147"
    print *, "    β = 0.693147 / 9.0e-4 ≈ 770.16"
    print *, "  Computed β = ", p%beta
    print *, ""
    
    ! 设置测试网格
    n = 1601  ! 高分辨率以显示细节
    allocate(x_vals(n), u_vals(n))
    
    ! 生成x坐标 [-1, 1]
    do i = 1, n
        x_vals(i) = -1.0_dp + 2.0_dp * (i-1) / (n-1)
    end do
    
    ! 计算初始条件
    print *, "=== Computing Initial Condition ==="
    do i = 1, n
        u_vals(i) = complex_wave_precise(x_vals(i))
    end do
    
    ! 保存结果
    call save_results(x_vals, u_vals, n, "precise_wave.dat")
    
    ! 关键点分析
    call analyze_key_points(x_vals, u_vals, n, p)
    
    ! β的意义分析
    call analyze_beta_meaning(p)
    
    deallocate(x_vals, u_vals)
    
    print *, ""
    print *, "=== Test completed ==="
    
contains
    
    subroutine save_results(x, u, n, filename)
        real(dp), intent(in) :: x(:), u(:)
        integer, intent(in) :: n
        character(len=*), intent(in) :: filename
        integer :: i
        
        open(20, file=filename, status='replace')
        write(20, *) "# Complex wave with precise parameters"
        write(20, *) "# a=0.5, z=-0.7, δ=0.005, α=10, β=log2/(36δ²)"
        write(20, *) "# x u"
        do i = 1, n
            write(20, '(2ES16.8)') x(i), u(i)
        end do
        close(20)
        print *, "  Saved to ", trim(filename)
    end subroutine
    
    subroutine analyze_key_points(x, u, n, p)
        real(dp), intent(in) :: x(:), u(:)
        integer, intent(in) :: n
        type(PreciseParams), intent(in) :: p
        
        integer :: i, idx
        real(dp) :: x_test, u_test, expected
        real(dp) :: g1, g2, g3, f1, f2, f3
        real(dp) :: arg1, arg2, arg3
        
        print *, ""
        print *, "=== Analysis at Key Points ==="
        
        ! 1. 区域1中心 (x = z = -0.7)
        x_test = -0.7_dp
        idx = minloc(abs(x - x_test), dim=1)
        u_test = u(idx)
        
        ! 理论计算
        g1 = exp(-p%beta * (x_test - (p%z - p%delta))**2)  ! G(z,β,z-δ)
        g2 = exp(-p%beta * (x_test - (p%z + p%delta))**2)  ! G(z,β,z+δ)
        g3 = exp(-p%beta * (x_test - p%z)**2)              ! G(z,β,z) = 1
        expected = (g1 + g2 + 4.0_dp * g3) / 6.0_dp
        
        print *, "1. Region 1 center (x = z = -0.7):"
        print *, "   G(z,β,z-δ) = exp(-βδ²) = ", g1
        print *, "   G(z,β,z+δ) = exp(-βδ²) = ", g2
        print *, "   G(z,β,z) = 1.0"
        print *, "   Expected u = (", g1, " + ", g2, " + 4*1.0) / 6 = ", expected
        print *, "   Actual u = ", u_test
        print *, "   Difference = ", abs(u_test - expected)
        print *, ""
        
        ! 2. 点 x = z ± δ
        x_test = p%z + p%delta  ! -0.695
        idx = minloc(abs(x - x_test), dim=1)
        u_test = u(idx)
        
        g1 = exp(-p%beta * (x_test - (p%z - p%delta))**2)  ! G(z+δ,β,z-δ)
        g2 = exp(-p%beta * (x_test - (p%z + p%delta))**2)  ! G(z+δ,β,z+δ) = 1
        g3 = exp(-p%beta * (x_test - p%z)**2)              ! G(z+δ,β,z)
        expected = (g1 + g2 + 4.0_dp * g3) / 6.0_dp
        
        print *, "2. Point x = z + δ = ", x_test, ":"
        print *, "   G(z+δ,β,z-δ) = exp(-β(2δ)²) = exp(-4βδ²) = ", g1
        print *, "   G(z+δ,β,z+δ) = 1.0"
        print *, "   G(z+δ,β,z) = exp(-βδ²) = ", g3
        print *, "   Expected u = ", expected
        print *, "   Actual u = ", u_test
        print *, ""
        
        ! 3. 区域4中心 (x = a = 0.5)
        x_test = 0.5_dp
        idx = minloc(abs(x - x_test), dim=1)
        u_test = u(idx)
        
        ! 理论计算
        arg1 = 1.0_dp - (p%alpha * (x_test - (p%a - p%delta)))**2
        arg2 = 1.0_dp - (p%alpha * (x_test - (p%a + p%delta)))**2
        arg3 = 1.0_dp - (p%alpha * (x_test - p%a))**2
        
        if (arg1 > 0.0_dp) then
            f1 = sqrt(arg1)
        else
            f1 = 0.0_dp
        end if
        
        if (arg2 > 0.0_dp) then
            f2 = sqrt(arg2)
        else
            f2 = 0.0_dp
        end if
        
        if (arg3 > 0.0_dp) then
            f3 = sqrt(arg3)  ! = 1.0
        else
            f3 = 0.0_dp
        end if
        
        expected = (f1 + f2 + 4.0_dp * f3) / 6.0_dp
        
        print *, "3. Region 4 center (x = a = 0.5):"
        print *, "   F(a,α,a-δ) = sqrt(1 - α²δ²) = ", f1
        print *, "   F(a,α,a+δ) = sqrt(1 - α²δ²) = ", f2
        print *, "   F(a,α,a) = 1.0"
        print *, "   Expected u = (", f1, " + ", f2, " + 4*1.0) / 6 = ", expected
        print *, "   Actual u = ", u_test
        print *, "   Difference = ", abs(u_test - expected)
        
    end subroutine
    
    subroutine analyze_beta_meaning(p)
        type(PreciseParams), intent(in) :: p
        
        real(dp) :: fwhm, sigma
        
        print *, ""
        print *, "=== Physical Meaning of β ==="
        print *, "For Gaussian: G(x) = exp(-β(x-z)²)"
        print *, ""
        
        ! 半高宽 (Full Width at Half Maximum)
        ! exp(-β*(FWHM/2)²) = 1/2
        ! -β*(FWHM/2)² = ln(1/2) = -ln(2)
        ! (FWHM/2)² = ln(2)/β
        ! FWHM = 2 * sqrt(ln(2)/β)
        fwhm = 2.0_dp * sqrt(log(2.0_dp) / p%beta)
        
        ! 标准差 σ (如果看作高斯分布：exp(-(x-z)²/(2σ²)))
        ! 比较：exp(-β(x-z)²) = exp(-(x-z)²/(2σ²))
        ! 所以 1/(2σ²) = β，σ = 1/sqrt(2β)
        sigma = 1.0_dp / sqrt(2.0_dp * p%beta)
        
        print *, "Gaussian properties:"
        print *, "  β = ", p%beta
        print *, "  Standard deviation σ = 1/√(2β) = ", sigma
        print *, "  Full Width at Half Maximum (FWHM) = 2√(ln2/β) = ", fwhm
        print *, "  In terms of δ: FWHM ≈ ", fwhm/p%delta, "δ"
        print *, ""
        
        ! 检查在 x = z ± δ 处的值
        print *, "At x = z ± δ:"
        print *, "  G(z±δ,β,z) = exp(-βδ²)"
        print *, "  Since β = ln2/(36δ²), βδ² = ln2/36 ≈ 0.01925"
        print *, "  So G(z±δ,β,z) = exp(-ln2/36) = 2^(-1/36) ≈ ", exp(-log(2.0_dp)/36.0_dp)
        print *, "  This is very close to 1.0 (about 0.981)"
        
    end subroutine
    
end program test_precise_parameters