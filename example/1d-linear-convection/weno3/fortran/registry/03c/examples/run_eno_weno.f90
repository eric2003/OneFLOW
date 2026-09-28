! examples/run_eno_weno.f90 (修复版)
program run_eno_weno
    ! 示例程序：ENO/WENO对比分析
    
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction, config_print
    use mesh_module, only: mesh_type
    use registry_module, only: registry_init, registry_cleanup, initialize_default_components
    use physics_solver_module, only: physics_solver, SOLVER_INITIALIZED, &
                                     SOLVER_COMPLETED, SOLVER_ERROR
    
    implicit none
    
    ! 求解器实例
    type(cfd_config) :: config_eno3, config_weno3, config_weno5
    type(mesh_type) :: mesh
    type(physics_solver) :: solver_eno3, solver_weno3, solver_weno5
    
    ! 保存结果的变量
    real(wp) :: time_eno3, time_weno3, time_weno5
    integer :: steps_eno3, steps_weno3, steps_weno5
    integer :: state_eno3, state_weno3, state_weno5
    character(len=100) :: error_eno3, error_weno3, error_weno5
    
    character(len=100) :: status_str
    logical :: all_success
    
    print *, "=========================================="
    print *, "  ENO/WENO对比分析 - Fortran版本"
    print *, "=========================================="
    print *, ""
    
    ! 步骤0: 系统初始化
    print *, "[步骤0] 初始化系统..."
    print *, "---------------------"
    
    ! 初始化注册系统
    call registry_init(verbose=.true.)
    
    ! 注册默认组件
    print *, "注册默认组件..."
    call initialize_default_components()
    print *, ""
    
    ! 步骤1: 创建网格
    print *, "[步骤1] 创建计算网格..."
    print *, "------------------------"
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=40)
    call mesh%print_info()
    print *, ""
    
    ! ========== ENO3 求解器 ==========
    print *, "[步骤2] 配置并运行ENO3求解器..."
    print *, "--------------------------------"
    
    ! 配置ENO3
    config_eno3%verbose = .true.
    config_eno3%ic_type = "step"
    config_eno3%wave_speed = 1.0_wp
    config_eno3%final_time = 0.625_wp
    config_eno3%dt = 0.0025_wp
    config_eno3%rk_order = 2
    config_eno3%boundary_type = "periodic"
    config_eno3%equation_type = "linear_advection"
    config_eno3%problem_type = "linear_advection"
    config_eno3%enable_physics = .true.
    config_eno3%domain_length = 2.0_wp
    
    call config_with_reconstruction(config_eno3, "eno", 3)
    
    print *, "ENO3配置信息:"
    call config_print(config_eno3)
    print *, ""
    
    ! 创建并运行ENO3求解器
    print *, "创建ENO3求解器实例..."
    solver_eno3 = physics_solver(config_eno3, mesh)
    
    print *, "初始化ENO3求解器..."
    call solver_eno3%initialize()
    
    print *, "运行ENO3求解器..."
    call solver_eno3%run_to_time(config_eno3%final_time)
    
    ! 立即保存ENO3结果
    state_eno3 = solver_eno3%get_state()
    time_eno3 = solver_eno3%current_time
    steps_eno3 = solver_eno3%current_step
    error_eno3 = solver_eno3%get_error()
    
    print *, "检查ENO3求解器状态..."
    if (state_eno3 == SOLVER_COMPLETED) then
        print *, "✓ ENO3求解器运行完成"
        print *, "  最终时间: ", time_eno3
        print *, "  总步数: ", steps_eno3
    else if (state_eno3 == SOLVER_ERROR) then
        print *, "✗ ENO3求解器运行失败"
        print *, "  错误信息: ", trim(error_eno3)
    else
        print *, "? ENO3求解器状态: ", state_eno3
    end if
    print *, ""
    
    ! ========== WENO3 求解器 ==========
    print *, "[步骤3] 配置并运行WENO3求解器..."
    print *, "---------------------------------"
    
    ! 配置WENO3
    config_weno3%verbose = .true.
    config_weno3%ic_type = "step"
    config_weno3%wave_speed = 1.0_wp
    config_weno3%final_time = 0.625_wp
    config_weno3%dt = 0.0025_wp
    config_weno3%rk_order = 2
    config_weno3%boundary_type = "periodic"
    config_weno3%equation_type = "linear_advection"
    config_weno3%problem_type = "linear_advection"
    config_weno3%enable_physics = .true.
    config_weno3%domain_length = 2.0_wp
    
    call config_with_reconstruction(config_weno3, "weno3", 3)
    
    print *, "WENO3配置信息:"
    call config_print(config_weno3)
    print *, ""
    
    ! 创建并运行WENO3求解器
    print *, "创建WENO3求解器实例..."
    solver_weno3 = physics_solver(config_weno3, mesh)
    
    print *, "初始化WENO3求解器..."
    call solver_weno3%initialize()
    
    print *, "运行WENO3求解器..."
    call solver_weno3%run_to_time(config_weno3%final_time)
    
    ! 立即保存WENO3结果
    state_weno3 = solver_weno3%get_state()
    time_weno3 = solver_weno3%current_time
    steps_weno3 = solver_weno3%current_step
    error_weno3 = solver_weno3%get_error()
    
    print *, "检查WENO3求解器状态..."
    if (state_weno3 == SOLVER_COMPLETED) then
        print *, "✓ WENO3求解器运行完成"
        print *, "  最终时间: ", time_weno3
        print *, "  总步数: ", steps_weno3
    else if (state_weno3 == SOLVER_ERROR) then
        print *, "✗ WENO3求解器运行失败"
        print *, "  错误信息: ", trim(error_weno3)
    else
        print *, "? WENO3求解器状态: ", state_weno3
    end if
    print *, ""
    
    ! ========== WENO5 求解器 ==========
    print *, "[步骤4] 配置并运行WENO5求解器..."
    print *, "---------------------------------"
    
    ! 配置WENO5
    config_weno5%verbose = .true.
    config_weno5%ic_type = "step"
    config_weno5%wave_speed = 1.0_wp
    config_weno5%final_time = 0.625_wp
    config_weno5%dt = 0.0025_wp
    config_weno5%rk_order = 2
    config_weno5%boundary_type = "periodic"
    config_weno5%equation_type = "linear_advection"
    config_weno5%problem_type = "linear_advection"
    config_weno5%enable_physics = .true.
    config_weno5%domain_length = 2.0_wp
    
    call config_with_reconstruction(config_weno5, "weno", 5)
    
    print *, "WENO5配置信息:"
    call config_print(config_weno5)
    print *, ""
    
    ! 创建并运行WENO5求解器
    print *, "创建WENO5求解器实例..."
    solver_weno5 = physics_solver(config_weno5, mesh)
    
    print *, "初始化WENO5求解器..."
    call solver_weno5%initialize()
    
    print *, "运行WENO5求解器..."
    call solver_weno5%run_to_time(config_weno5%final_time)
    
    ! 立即保存WENO5结果
    state_weno5 = solver_weno5%get_state()
    time_weno5 = solver_weno5%current_time
    steps_weno5 = solver_weno5%current_step
    error_weno5 = solver_weno5%get_error()
    
    print *, "检查WENO5求解器状态..."
    if (state_weno5 == SOLVER_COMPLETED) then
        print *, "✓ WENO5求解器运行完成"
        print *, "  最终时间: ", time_weno5
        print *, "  总步数: ", steps_weno5
    else if (state_weno5 == SOLVER_ERROR) then
        print *, "✗ WENO5求解器运行失败"
        print *, "  错误信息: ", trim(error_weno5)
    else
        print *, "? WENO5求解器状态: ", state_weno5
    end if
    print *, ""
    
    ! ========== 清理求解器 ==========
    print *, "[步骤5] 清理求解器..."
    print *, "----------------------"
    
    call solver_eno3%cleanup()
    call solver_weno3%cleanup()
    call solver_weno5%cleanup()
    
    ! ========== 清理注册系统 ==========
    print *, "[步骤6] 清理注册系统..."
    print *, "----------------------"
    call registry_cleanup()
    
    ! ========== 结果汇总 ==========
    print *, ""
    print *, "=========================================="
    print *, "            结果汇总"
    print *, "=========================================="
    print *, ""
    
    print *, "求解器状态对比:"
    print *, "---------------"
    
    ! ENO3状态
    if (state_eno3 == SOLVER_COMPLETED) then
        status_str = "完成 ✓"
    else if (state_eno3 == SOLVER_INITIALIZED) then
        status_str = "已初始化"
    else if (state_eno3 == SOLVER_ERROR) then
        status_str = "错误 ✗"
    else
        write(status_str, '(A, I3)') "未知状态 ", state_eno3
    end if
    print *, "ENO3:   ", trim(status_str)
    print *, "        最终时间: ", time_eno3
    print *, "        总步数:   ", steps_eno3
    if (len_trim(error_eno3) > 0) then
        print *, "        错误信息: ", trim(error_eno3)
    end if
    print *, ""
    
    ! WENO3状态
    if (state_weno3 == SOLVER_COMPLETED) then
        status_str = "完成 ✓"
    else if (state_weno3 == SOLVER_INITIALIZED) then
        status_str = "已初始化"
    else if (state_weno3 == SOLVER_ERROR) then
        status_str = "错误 ✗"
    else
        write(status_str, '(A, I3)') "未知状态 ", state_weno3
    end if
    print *, "WENO3:  ", trim(status_str)
    print *, "        最终时间: ", time_weno3
    print *, "        总步数:   ", steps_weno3
    if (len_trim(error_weno3) > 0) then
        print *, "        错误信息: ", trim(error_weno3)
    end if
    print *, ""
    
    ! WENO5状态
    if (state_weno5 == SOLVER_COMPLETED) then
        status_str = "完成 ✓"
    else if (state_weno5 == SOLVER_INITIALIZED) then
        status_str = "已初始化"
    else if (state_weno5 == SOLVER_ERROR) then
        status_str = "错误 ✗"
    else
        write(status_str, '(A, I3)') "未知状态 ", state_weno5
    end if
    print *, "WENO5:  ", trim(status_str)
    print *, "        最终时间: ", time_weno5
    print *, "        总步数:   ", steps_weno5
    if (len_trim(error_weno5) > 0) then
        print *, "        错误信息: ", trim(error_weno5)
    end if
    print *, ""
    
    ! ========== 最终判断 ==========
    all_success = (state_eno3 == SOLVER_COMPLETED) .and. &
                  (state_weno3 == SOLVER_COMPLETED) .and. &
                  (state_weno5 == SOLVER_COMPLETED)
    
    if (all_success) then
        print *, "✓ 所有求解器成功运行！"
        print *, ""
        print *, "计算参数总结:"
        print *, "--------------"
        print *, "网格单元数:     ", mesh%ncells
        print *, "时间步长:       ", config_eno3%dt
        print *, "最终时间:       ", config_eno3%final_time
        print *, "波速:           ", config_eno3%wave_speed
        print *, "边界条件:       ", trim(config_eno3%boundary_type)
        print *, "初始条件:       ", trim(config_eno3%ic_type)
        print *, ""
        print *, "重构格式对比:"
        print *, "--------------"
        print *, "ENO3:   3阶本质无振荡"
        print *, "WENO3:  3阶加权本质无振荡"
        print *, "WENO5:  5阶加权本质无振荡"
        print *, ""
        print *, "性能对比:"
        print *, "----------"
        print *, "ENO3:  ", steps_eno3, " 步，最终时间 ", time_eno3
        print *, "WENO3: ", steps_weno3, " 步，最终时间 ", time_weno3
        print *, "WENO5: ", steps_weno5, " 步，最终时间 ", time_weno5
        print *, ""
        print *, "下一步开发计划:"
        print *, "1. 实现真实的ENO/WENO重构算法"
        print *, "2. 实现Rusanov通量计算"
        print *, "3. 添加RK时间积分器"
        print *, "4. 实现结果输出和可视化"
    else
        print *, "✗ 有求解器运行失败"
        print *, ""
        print *, "故障排除:"
        print *, "----------"
        
        if (state_eno3 /= SOLVER_COMPLETED) then
            print *, "• ENO3失败: 状态=", state_eno3
            if (len_trim(error_eno3) > 0) then
                print *, "  错误信息: ", trim(error_eno3)
            end if
        end if
        
        if (state_weno3 /= SOLVER_COMPLETED) then
            print *, "• WENO3失败: 状态=", state_weno3
            if (len_trim(error_weno3) > 0) then
                print *, "  错误信息: ", trim(error_weno3)
            end if
        end if
        
        if (state_weno5 /= SOLVER_COMPLETED) then
            print *, "• WENO5失败: 状态=", state_weno5
            if (len_trim(error_weno5) > 0) then
                print *, "  错误信息: ", trim(error_weno5)
            end if
        end if
        
        print *, ""
        print *, "可能的原因:"
        print *, "1. 配置参数不正确"
        print *, "2. 内存分配失败"
        print *, "3. 数值计算不稳定"
    end if
    
    print *, ""
    print *, "=========================================="
    print *, "              分析完成"
    print *, "=========================================="
    
end program run_eno_weno