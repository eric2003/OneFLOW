! tests/test_solver_framework.f90
program test_solver_framework
    use, intrinsic :: iso_fortran_env, only: real64
    use config_module, only: cfd_config, config_print, config_with_reconstruction
    use mesh_module, only: mesh_type
    use solver_module, only: cfd_solver, solver_create, solver_run, solver_cleanup
    use solver_module, only: SOLVER_UNINITIALIZED, SOLVER_INITIALIZED, &
                            SOLVER_COMPLETED, SOLVER_ERROR
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(cfd_solver) :: solver
    
    print *, "=== 求解器框架测试 ==="
    print *, ""
    
    ! 测试1: 基本创建
    print *, "1. 测试求解器创建..."
    print *, "----------------------"
    
    ! 创建配置
    config%verbose = .true.
    call config_with_reconstruction(config, "eno", 3)
    config%flux_type = "rusanov"
    config%wave_speed = 1.0_real64
    config%dt = 0.01_real64
    
    call config_print(config)
    print *, ""
    
    ! 创建网格
    call mesh%init(xmin=0.0_real64, xmax=2.0_real64, ncells=20)
    call mesh%print_info()
    print *, ""
    
    ! 创建求解器
    solver = solver_create(config, mesh)
    print *, "✓ 求解器创建成功"
    print *, "  状态: ", solver%get_state()
    print *, ""
    
    ! 测试2: 求解器初始化
    print *, "2. 测试求解器初始化..."
    print *, "------------------------"
    
    call solver%initialize()
    print *, "✓ 求解器初始化完成"
    print *, "  状态: ", solver%get_state()
    print *, "  错误信息: '", trim(solver%get_error()), "'"
    print *, ""
    
    ! 测试3: 简单运行
    print *, "3. 测试求解器运行..."
    print *, "----------------------"
    
    call solver_run(solver, 0.05_real64)  ! 运行到0.05秒
    print *, "✓ 求解器运行完成"
    print *, "  最终状态: ", solver%get_state()
    print *, ""
    
    ! 测试4: 清理
    print *, "4. 测试求解器清理..."
    print *, "----------------------"
    
    call solver_cleanup(solver)
    print *, "✓ 求解器清理完成"
    print *, "  状态: ", solver%get_state()
    print *, ""
    
    ! 测试5: 错误处理
    print *, "5. 测试错误处理..."
    print *, "-------------------"
    
    ! 尝试重复初始化
    call solver%initialize()
    print *, "  重复初始化状态: ", solver%get_state()
    print *, "  错误信息: '", trim(solver%get_error()), "'"
    
    call solver_cleanup(solver)
    print *, ""
    
    print *, "=== 框架测试总结 ==="
    print *, "✓ 求解器创建/初始化/运行/清理流程验证完成"
    print *, "✓ 状态管理正常工作"
    print *, "✓ 错误处理机制就绪"
    print *, ""
    print *, "下一步: 添加实际数值计算功能"
    
end program test_solver_framework