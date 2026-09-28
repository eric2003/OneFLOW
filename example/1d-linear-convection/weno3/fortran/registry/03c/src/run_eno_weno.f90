! src/run_eno_weno.f90
program run_eno_weno
    ! 使用所有需要的模块
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction, config_print
    use mesh_module, only: mesh_type
    use solver_base_module, only: SOLVER_INITIALIZED, SOLVER_COMPLETED
    
    use physics_solver_module, only: physics_solver
    use registry_module, only: registry_init, registry_cleanup
    
    implicit none
    
    type(cfd_config) :: config_eno3, config_weno3, config_weno5
    type(mesh_type) :: mesh
    type(physics_solver) :: solver_eno3, solver_weno3, solver_weno5
    
    character(len=100) :: output_file
    integer :: status
    
    print *, "=== ENO/WENO对比分析 (Fortran版本) ==="
    print *, ""
    
    ! 初始化注册系统
    call registry_init(verbose=.true.)
    print *, ""
    
    ! 步骤1: 初始化网格
    print *, "步骤1: 初始化网格..."
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=40)
    call mesh%print_info()
    print *, ""
    
    ! 步骤2: 配置并运行ENO3求解器
    print *, "步骤2: 运行ENO3求解器..."
    print *, "-------------------------"
    
    ! 创建ENO3配置
    config_eno3%verbose = .true.
    call config_with_reconstruction(config_eno3, "eno", 3)
    config_eno3%dt = 0.0025_wp
    config_eno3%rk_order = 2
    config_eno3%wave_speed = 1.0_wp
    config_eno3%final_time = 0.625_wp
    config_eno3%ic_type = "step"
    
    print *, "ENO3配置:"
    call config_print(config_eno3)
    
    ! 创建ENO3求解器
    print *, "创建ENO3求解器..."
    solver_eno3 = physics_solver(config_eno3, mesh)
    call solver_eno3%print_info()
    print *, ""
    
    ! 初始化并运行ENO3求解器
    print *, "运行ENO3求解器..."
    call solver_eno3%initialize()
    call solver_eno3%run_to_time(config_eno3%final_time)
    
    if (solver_eno3%get_state() == SOLVER_COMPLETED) then
        print *, "✓ ENO3求解器运行完成"
        print *, "  最终时间: ", solver_eno3%current_time
        print *, "  总步数: ", solver_eno3%current_step
    else
        print *, "✗ ENO3求解器运行失败"
        print *, "  错误信息: ", trim(solver_eno3%get_error())
    end if
    print *, ""
    
    ! 步骤3: 配置并运行WENO3求解器
    print *, "步骤3: 运行WENO3求解器..."
    print *, "--------------------------"
    
    ! 创建WENO3配置
    config_weno3%verbose = .true.
    call config_with_reconstruction(config_weno3, "weno3", 3)
    config_weno3%dt = 0.0025_wp
    config_weno3%rk_order = 2
    config_weno3%wave_speed = 1.0_wp
    config_weno3%final_time = 0.625_wp
    config_weno3%ic_type = "step"
    
    print *, "WENO3配置:"
    call config_print(config_weno3)
    
    ! 创建WENO3求解器
    print *, "创建WENO3求解器..."
    solver_weno3 = physics_solver(config_weno3, mesh)
    call solver_weno3%print_info()
    print *, ""
    
    ! 初始化并运行WENO3求解器
    print *, "运行WENO3求解器..."
    call solver_weno3%initialize()
    call solver_weno3%run_to_time(config_weno3%final_time)
    
    if (solver_weno3%get_state() == SOLVER_COMPLETED) then
        print *, "✓ WENO3求解器运行完成"
        print *, "  最终时间: ", solver_weno3%current_time
        print *, "  总步数: ", solver_weno3%current_step
    else
        print *, "✗ WENO3求解器运行失败"
        print *, "  错误信息: ", trim(solver_weno3%get_error())
    end if
    print *, ""
    
    ! 步骤4: 配置并运行WENO5求解器
    print *, "步骤4: 运行WENO5求解器..."
    print *, "--------------------------"
    
    ! 创建WENO5配置
    config_weno5%verbose = .true.
    call config_with_reconstruction(config_weno5, "weno", 5)
    config_weno5%dt = 0.0025_wp
    config_weno5%rk_order = 2
    config_weno5%wave_speed = 1.0_wp
    config_weno5%final_time = 0.625_wp
    config_weno5%ic_type = "step"
    
    print *, "WENO5配置:"
    call config_print(config_weno5)
    
    ! 创建WENO5求解器
    print *, "创建WENO5求解器..."
    solver_weno5 = physics_solver(config_weno5, mesh)
    call solver_weno5%print_info()
    print *, ""
    
    ! 初始化并运行WENO5求解器
    print *, "运行WENO5求解器..."
    call solver_weno5%initialize()
    call solver_weno5%run_to_time(config_weno5%final_time)
    
    if (solver_weno5%get_state() == SOLVER_COMPLETED) then
        print *, "✓ WENO5求解器运行完成"
        print *, "  最终时间: ", solver_weno5%current_time
        print *, "  总步数: ", solver_weno5%current_step
    else
        print *, "✗ WENO5求解器运行失败"
        print *, "  错误信息: ", trim(solver_weno5%get_error())
    end if
    print *, ""
    
    ! 清理注册系统
    call registry_cleanup()
    
    ! 清理求解器
    call solver_eno3%cleanup()
    call solver_weno3%cleanup()
    call solver_weno5%cleanup()
    
    print *, "=== 分析完成 ==="
    print *, "总结:"
    print *, "  ENO3:  时间=", solver_eno3%current_time, " 步数=", solver_eno3%current_step
    print *, "  WENO3: 时间=", solver_weno3%current_time, " 步数=", solver_weno3%current_step
    print *, "  WENO5: 时间=", solver_weno5%current_time, " 步数=", solver_weno5%current_step
    
    if (solver_eno3%get_state() == SOLVER_COMPLETED .and. &
        solver_weno3%get_state() == SOLVER_COMPLETED .and. &
        solver_weno5%get_state() == SOLVER_COMPLETED) then
        print *, ""
        print *, "✓ 所有求解器成功运行！"
        print *, "下一步: 添加边界条件模块"
    else
        print *, ""
        print *, "✗ 有求解器运行失败"
        print *, "需要先调试现有代码"
    end if
    
end program run_eno_weno