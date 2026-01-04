! tests/test_physics_minimal.f90
program test_physics_minimal
    use precision_module, only: wp, ip
    use linear_convection_equation, only: linear_convection_eq, create_linear_convection_eq
    use linear_convection_problem, only: linear_convection_prob, create_linear_convection_prob
    
    implicit none
    
    type(linear_convection_eq) :: eq
    type(linear_convection_prob) :: prob
    real(wp) :: u, f, a
    real(wp), allocatable :: x(:), u_ic(:), u_exact(:)
    integer :: i, nx = 10
    
    print *, "=== 最小物理模块测试 ==="
    print *, ""
    
    ! 测试1: 方程功能
    print *, "1. 测试方程功能..."
    print *, "-------------------"
    
    eq = create_linear_convection_eq(wave_speed=2.0_wp)
    print *, "方程: ", eq%name
    print *, "波速: ", eq%wave_speed
    
    u = 1.5_wp
    f = eq%flux(u)
    a = eq%speed()
    
    print *, "u = ", u
    print *, "F(u) = ", f, " (期望: 3.0)"
    print *, "波速 a = ", a, " (期望: 2.0)"
    
    if (abs(f - 3.0_wp) < 1e-10_wp .and. abs(a - 2.0_wp) < 1e-10_wp) then
        print *, "✓ 方程功能正常"
    else
        print *, "✗ 方程功能异常"
    end if
    print *, ""
    
    ! 测试2: 问题功能
    print *, "2. 测试问题功能..."
    print *, "-------------------"
    
    prob = create_linear_convection_prob(ic_type="step", domain_length=2.0_wp)
    print *, "问题: ", prob%name
    print *, "IC类型: ", trim(prob%ic_type)
    print *, "域长度: ", prob%domain_length
    
    allocate(x(nx), u_ic(nx), u_exact(nx))
    do i = 1, nx
        x(i) = 0.0_wp + (i-1) * 0.2_wp
    end do
    
    ! 测试初始条件
    call prob%initial_condition(x, u_ic)
    print *, "初始条件范围: ", minval(u_ic), " 到 ", maxval(u_ic)
    
    ! 测试精确解
    u_exact = prob%exact_solution(x, 0.0_wp)
    print *, "t=0时精确解范围: ", minval(u_exact), " 到 ", maxval(u_exact)
    
    ! 检查阶跃函数
    if (abs(u_ic(1) - 1.0_wp) < 1e-10_wp .and. &
        abs(u_ic(6) - 2.0_wp) < 1e-10_wp) then
        print *, "✓ 阶跃初始条件正确"
    else
        print *, "✗ 阶跃初始条件错误"
    end if
    
    ! 检查精确解与初始条件一致
    if (maxval(abs(u_ic - u_exact)) < 1e-10_wp) then
        print *, "✓ t=0时精确解与初始条件一致"
    else
        print *, "✗ 精确解计算错误"
    end if
    print *, ""
    
    ! 测试3: 边界条件接口
    print *, "3. 测试边界条件接口..."
    print *, "----------------------"
    call prob%boundary_condition(u_ic, 0.0_wp)
    print *, "✓ 边界条件接口正常"
    print *, ""
    
    deallocate(x, u_ic, u_exact)
    
    print *, "=== 物理模块最小测试完成 ==="
    print *, "下一步: 将物理模块集成到现有系统中"
    
end program test_physics_minimal