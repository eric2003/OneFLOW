! tests/test_physics_solver_simple.f90
program test_physics_solver_simple
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction
    use mesh_module, only: mesh_type
    use physics_solver_module, only: physics_solver, SOLVER_COMPLETED
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(physics_solver) :: solver
    real(wp) :: final_time, final_step
    integer :: state
    
    print *, "========================================="
    print *, "  简单物理求解器测试"
    print *, "========================================="
    print *, ""
    
    ! 步骤1: 配置
    print *, "[步骤1] 配置求解器..."
    print *, "---------------------"
    
    config%verbose = .true.
    call config_with_reconstruction(config, "eno", 3)
    config%dt = 0.01_wp
    config%final_time = 0.1_wp
    config%wave_speed = 1.0_wp
    config%ic_type = "step"
    config%boundary_type = "periodic"
    config%equation_type = "linear_advection"
    config%problem_type = "linear_advection"
    config%enable_physics = .true.
    config%domain_length = 1.0_wp
    
    print *, "配置参数:"
    print *, "  重构格式: ", trim(config%recon_scheme)
    print *, "  时间步长: ", config%dt
    print *, "  最终时间: ", config%final_time
    print *, "  波速: ", config%wave_speed
    print *, ""
    
    ! 步骤2: 创建网格
    print *, "[步骤2] 创建网格..."
    print *, "-------------------"
    
    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=10)
    
    print *, "网格信息:"
    print *, "  单元数: ", mesh%ncells
    print *, "  节点数: ", mesh%nnodes
    print *, "  网格间距: ", mesh%dx
    print *, ""
    
    ! 步骤3: 创建求解器
    print *, "[步骤3] 创建求解器..."
    print *, "---------------------"
    
    solver = physics_solver(config, mesh)
    
    print *, "求解器创建成功"
    print *, "  初始状态: ", solver%get_state()
    print *, ""
    
    ! 步骤4: 初始化
    print *, "[步骤4] 初始化求解器..."
    print *, "-----------------------"
    
    call solver%initialize()
    
    state = solver%get_state()
    print *, "初始化完成"
    print *, "  状态: ", state
    print *, "  当前时间: ", solver%current_time
    print *, "  当前步数: ", solver%current_step
    print *, ""
    
    ! 步骤5: 运行求解器
    print *, "[步骤5] 运行求解器..."
    print *, "---------------------"
    
    call solver%run_to_time(config%final_time)
    
    state = solver%get_state()
    print *, "运行完成"
    print *, "  状态: ", state
    print *, "  最终时间: ", solver%current_time
    print *, "  总步数: ", solver%current_step
    print *, ""
    
    ! 步骤6: 保存结果
    print *, "[步骤6] 保存结果..."
    print *, "-------------------"
    
    final_time = solver%current_time
    final_step = real(solver%current_step, wp)
    state = solver%get_state()
    
    print *, "保存的结果:"
    print *, "  状态: ", state
    print *, "  时间: ", final_time
    print *, "  步数: ", final_step
    print *, ""
    
    ! 步骤7: 清理求解器
    print *, "[步骤7] 清理求解器..."
    print *, "---------------------"
    
    call solver%cleanup()
    
    print *, "清理后状态:"
    print *, "  状态: ", solver%get_state()
    print *, "  时间: ", solver%current_time
    print *, "  步数: ", solver%current_step
    print *, ""
    
    ! 步骤8: 验证结果
    print *, "[步骤8] 验证结果..."
    print *, "-------------------"
    
    print *, "验证标准:"
    print *, "  1. 运行后状态应为 COMPLETED (", SOLVER_COMPLETED, ")"
    print *, "  2. 最终时间应接近 ", config%final_time
    print *, "  3. 步数应大于 0"
    print *, ""
    
    if (state == SOLVER_COMPLETED) then
        print *, "✓ 状态验证通过: COMPLETED"
    else
        print *, "✗ 状态验证失败: 期望 ", SOLVER_COMPLETED, ", 实际 ", state
    end if
    
    if (abs(final_time - config%final_time) < 1e-5_wp) then
        print *, "✓ 时间验证通过: ", final_time, " ≈ ", config%final_time
    else
        print *, "✗ 时间验证失败: ", final_time, " ≠ ", config%final_time
    end if
    
    if (final_step > 0) then
        print *, "✓ 步数验证通过: ", final_step, " > 0"
    else
        print *, "✗ 步数验证失败: ", final_step, " ≤ 0"
    end if
    
    print *, ""
    
    ! 最终判断
    if (state == SOLVER_COMPLETED .and. &
        abs(final_time - config%final_time) < 1e-5_wp .and. &
        final_step > 0) then
        print *, "========================================="
        print *, "  所有测试通过！ ✓"
        print *, "========================================="
    else
        print *, "========================================="
        print *, "  测试失败 ✗"
        print *, "========================================="
    end if
    
end program test_physics_solver_simple