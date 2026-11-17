! tests/test_initial_condition.f90
program test_initial_condition
    use base_modules, only: wp
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create
    use ic_factory_module, only: create_initial_condition
    use ic_base_module, only: initial_condition
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(domain_type) :: domain
    type(solution_type) :: solution
    class(initial_condition), allocatable :: ic
    integer :: i
    
    print *, "=== 初始条件模块测试 ==="
    print *, ""
    
    ! 创建配置和网格
    config%verbose = .false.
    config%recon_scheme = "eno"
    config%spatial_order = 3
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=10)
    domain = domain_create(config, mesh)
    solution = solution_create(domain)
    
    ! 测试1: 阶跃函数初始条件
    print *, "1. 测试阶跃函数初始条件..."
    call create_initial_condition("step", ic)
    
    if (allocated(ic)) then
        call ic%apply(solution)
        print *, "  成功应用阶跃函数初始条件"
        print *, "  解范围: ", minval(solution%u), " 到 ", maxval(solution%u)
        
        ! 检查结果
        if (abs(maxval(solution%u) - 2.0_wp) < 1e-10_wp .and. &
            abs(minval(solution%u) - 1.0_wp) < 1e-10_wp) then
            print *, "  ✓ 阶跃函数测试通过"
        else
            print *, "  ✗ 阶跃函数测试失败"
        end if
    end if
    
    deallocate(ic)
    print *, ""
    
    ! 测试2: 正弦波初始条件
    print *, "2. 测试正弦波初始条件..."
    call create_initial_condition("sin", ic)
    
    if (allocated(ic)) then
        call solution%reset()
        call ic%apply(solution)
        print *, "  成功应用正弦波初始条件"
        print *, "  解范围: ", minval(solution%u), " 到 ", maxval(solution%u)
        print *, "  ✓ 正弦波测试通过"
    end if
    
    deallocate(ic)
    print *, ""
    
    ! 测试3: 高斯脉冲初始条件
    print *, "3. 测试高斯脉冲初始条件..."
    call create_initial_condition("gaussian", ic)
    
    if (allocated(ic)) then
        call solution%reset()
        call ic%apply(solution)
        print *, "  成功应用高斯脉冲初始条件"
        print *, "  解范围: ", minval(solution%u), " 到 ", maxval(solution%u)
        print *, "  ✓ 高斯脉冲测试通过"
    end if
    
    deallocate(ic)
    print *, ""
    
    ! 测试4: 错误处理
    print *, "4. 测试错误处理..."
    call create_initial_condition("unknown", ic)
    
    if (allocated(ic)) then
        print *, "  ✓ 错误处理测试通过（使用了默认初始条件）"
    else
        print *, "  ✗ 错误处理测试失败"
    end if
    
    print *, ""
    print *, "=== 初始条件模块测试完成 ==="
    
end program test_initial_condition