! examples/run_eno_weno_integrated.f90 (完整修复版)
program run_eno_weno_integrated
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction, config_print
    use mesh_module, only: mesh_type
    use solver_integrated_module, only: integrated_solver, SOLVER_INITIALIZED, SOLVER_COMPLETED
    
    implicit none
    
    type(cfd_config) :: config_eno3, config_weno3, config_weno5
    type(mesh_type) :: mesh
    type(integrated_solver) :: solver_eno3, solver_weno3, solver_weno5
    
    ! 结果变量
    real(wp) :: time_eno3, time_weno3, time_weno5
    integer :: steps_eno3, steps_weno3, steps_weno5
    integer :: state_eno3, state_weno3, state_weno5
    logical :: all_success
    
    ! 定义常量（避免导入冲突）
    integer, parameter :: SOLVER_READY = 0
    integer, parameter :: SOLVER_RUNNING = 2
    integer, parameter :: SOLVER_ERROR = -1
    
    print *, "=========================================="
    print *, "ENO/WENO 对比分析 (集成版本)"
    print *, "=========================================="
    print *, ""
    
    ! 步骤1: 创建网格
    print *, "[STEP 1] 创建计算网格..."
    print *, "-----------------------------------"
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=40)
    call mesh%print_info()
    print *, ""
    
    ! ========== ENO3 求解器 ==========
    print *, "[STEP 2] 配置和运行 ENO3 求解器..."
    print *, "-----------------------------------"
    
    ! 初始化配置
    config_eno3%verbose = .true.
    config_eno3%ic_type = "step"
    config_eno3%wave_speed = 1.0_wp
    config_eno3%final_time = 0.625_wp
    config_eno3%dt = 0.0025_wp
    config_eno3%boundary_type = "periodic"
    config_eno3%equation_type = "linear_advection"
    config_eno3%problem_type = "linear_advection"
    config_eno3%domain_length = 2.0_wp
    config_eno3%enable_physics = .true.
    
    call config_with_reconstruction(config_eno3, "eno", 3)
    call config_print(config_eno3)
    print *, ""
    
    ! 创建和运行求解器
    solver_eno3%config = config_eno3
    solver_eno3%mesh = mesh
    
    ! 控制数据模式 (当前使用简单数据)
    call solver_eno3%enable_real_data(.false.)  ! 设置为false使用简单数据
    
    call solver_eno3%initialize()
    call solver_eno3%run_to_time(config_eno3%final_time)
    
    ! 获取结果
    time_eno3 = solver_eno3%current_time
    steps_eno3 = solver_eno3%current_step
    state_eno3 = solver_eno3%get_state()
    
    print *, "ENO3 完成:"
    print *, "  最终时间: ", time_eno3
    print *, "  总步数: ", steps_eno3
    print *, "  状态: ", state_eno3
    print *, ""
    
    ! ========== WENO3 求解器 ==========
    print *, "[STEP 3] 配置和运行 WENO3 求解器..."
    print *, "-----------------------------------"
    
    ! 配置 (复制ENO3配置，只改重构格式)
    config_weno3 = config_eno3
    call config_with_reconstruction(config_weno3, "weno3", 3)
    
    ! 运行
    solver_weno3%config = config_weno3
    solver_weno3%mesh = mesh
    call solver_weno3%enable_real_data(.false.)
    
    call solver_weno3%initialize()
    call solver_weno3%run_to_time(config_weno3%final_time)
    
    time_weno3 = solver_weno3%current_time
    steps_weno3 = solver_weno3%current_step
    state_weno3 = solver_weno3%get_state()
    
    print *, "WENO3 完成:"
    print *, "  最终时间: ", time_weno3
    print *, "  总步数: ", steps_weno3
    print *, "  状态: ", state_weno3
    print *, ""
    
    ! ========== WENO5 求解器 ==========
    print *, "[STEP 4] 配置和运行 WENO5 求解器..."
    print *, "-----------------------------------"
    
    config_weno5 = config_eno3
    call config_with_reconstruction(config_weno5, "weno", 5)
    
    solver_weno5%config = config_weno5
    solver_weno5%mesh = mesh
    call solver_weno5%enable_real_data(.false.)
    
    call solver_weno5%initialize()
    call solver_weno5%run_to_time(config_weno5%final_time)
    
    time_weno5 = solver_weno5%current_time
    steps_weno5 = solver_weno5%current_step
    state_weno5 = solver_weno5%get_state()
    
    print *, "WENO5 完成:"
    print *, "  最终时间: ", time_weno5
    print *, "  总步数: ", steps_weno5
    print *, "  状态: ", state_weno5
    print *, ""
    
    ! ========== 结果汇总 ==========
    print *, "=========================================="
    print *, "          结果汇总"
    print *, "=========================================="
    print *, ""
    
    all_success = (state_eno3 == SOLVER_COMPLETED) .and. &
                  (state_weno3 == SOLVER_COMPLETED) .and. &
                  (state_weno5 == SOLVER_COMPLETED)
    
    if (all_success) then
        print *, "✓ 所有求解器成功完成!"
        print *, ""
        print *, "性能对比:"
        print *, "  ENO3:  ", steps_eno3, " 步"
        print *, "  WENO3: ", steps_weno3, " 步"
        print *, "  WENO5: ", steps_weno5, " 步"
        print *, ""
        print *, "网格信息:"
        print *, "  单元数: ", mesh%ncells
        print *, "  域长度: ", mesh%L
        print *, "  网格间距: ", mesh%dx
        print *, ""
        print *, "数据模式: 简单数据 (一阶迎风格式)"
        print *, "切换真实计算: 调用 solver%enable_real_data(.true.)"
    else
        print *, "✗ 部分求解器失败"
        print *, "  ENO3状态: ", state_eno3, " (期望: ", SOLVER_COMPLETED, ")"
        print *, "  WENO3状态: ", state_weno3, " (期望: ", SOLVER_COMPLETED, ")"
        print *, "  WENO5状态: ", state_weno5, " (期望: ", SOLVER_COMPLETED, ")"
    end if
    
    print *, ""
    print *, "=========================================="
    print *, "          分析完成"
    print *, "=========================================="
    
    ! 清理
    call solver_eno3%cleanup()
    call solver_weno3%cleanup()
    call solver_weno5%cleanup()
    
    ! 等待用户输入
    print *, ""
    print *, "按 ENTER 键退出..."
    read(*,*)
    
end program run_eno_weno_integrated