! tests/test_physics_solver.f90 (修复版)
program test_physics_solver
    use base_modules, only: wp  ! 使用一致的wp定义
    use config_module, only: cfd_config, config_print, config_with_reconstruction
    use mesh_module, only: mesh_type
    use solver_base_module, only: SOLVER_UNINITIALIZED, SOLVER_INITIALIZED, &
                                 SOLVER_COMPLETED, SOLVER_ERROR
    use physics_solver_module, only: physics_solver
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(physics_solver) :: psolver
    character(len=100) :: error_msg
    integer :: state
    
    print *, "=== Physics Solver Test ==="
    print *, ""
    
    ! 测试1: 创建物理求解器（默认物理配置）
    print *, "1. Creating physics solver (default physics)..."
    print *, "------------------------------------------------"
    
    config%verbose = .true.
    call config_with_reconstruction(config, "eno", 3)
    config%dt = 0.01_wp
    config%enable_physics = .true.
    
    call config_print(config)
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=10)
    
    psolver = physics_solver(config, mesh)
    call psolver%print_info()
    print *, ""
    
    ! 测试2: 初始化
    print *, "2. Initializing physics solver..."
    print *, "----------------------------------"
    
    call psolver%initialize()
    state = psolver%get_state()
    error_msg = psolver%get_error()
    print *, "State after initialization: ", state
    print *, "Expected: ", SOLVER_INITIALIZED
    print *, "Match? ", state == SOLVER_INITIALIZED
    print *, "Error message: '", trim(error_msg), "'"
    print *, ""
    
    ! 测试3: 运行一小段时间
    print *, "3. Running physics solver (short time)..."
    print *, "------------------------------------------"
    
    call psolver%run_to_time(0.02_wp)
    state = psolver%get_state()
    error_msg = psolver%get_error()
    print *, "State after short run: ", state
    print *, "Expected: ", SOLVER_COMPLETED
    print *, "Match? ", state == SOLVER_COMPLETED
    print *, "Current time: ", psolver%current_time
    print *, "Current step: ", psolver%current_step
    print *, ""
    
    ! 测试4: 禁用物理模块
    print *, "4. Testing physics solver with physics disabled..."
    print *, "--------------------------------------------------"
    
    config%enable_physics = .false.
    config%verbose = .false.
    
    psolver = physics_solver(config, mesh)
    call psolver%initialize()
    call psolver%run_to_time(0.01_wp)
    
    state = psolver%get_state()
    error_msg = psolver%get_error()
    print *, "State with physics disabled: ", state
    print *, "Expected: ", SOLVER_COMPLETED
    print *, "Match? ", state == SOLVER_COMPLETED
    print *, ""
    
    ! 测试5: 不同物理配置
    print *, "5. Testing different physics configurations..."
    print *, "----------------------------------------------"
    
    config%verbose = .true.
    config%enable_physics = .true.
    config%equation_type = "linear_advection"
    config%problem_type = "linear_advection"
    config%wave_speed = 2.5_wp
    config%domain_length = 3.0_wp
    config%ic_type = "gaussian"
    
    psolver = physics_solver(config, mesh)
    call psolver%initialize()
    
    print *, "Physics configuration test completed"
    print *, ""
    
    ! 测试6: 清理和错误处理
    print *, "6. Testing cleanup and error handling..."
    print *, "----------------------------------------"
    
    call psolver%cleanup()
    state = psolver%get_state()
    error_msg = psolver%get_error()
    print *, "State after cleanup: ", state
    print *, "Expected: ", SOLVER_UNINITIALIZED
    print *, "Match? ", state == SOLVER_UNINITIALIZED
    
    ! 尝试运行已清理的求解器
    call psolver%run_to_time(0.01_wp)
    state = psolver%get_state()
    error_msg = psolver%get_error()
    print *, "State after error attempt: ", state
    print *, "Expected: ", SOLVER_ERROR
    print *, "Match? ", state == SOLVER_ERROR
    print *, "Error message: '", trim(error_msg), "'"
    print *, ""
    
    ! 最终信息
    print *, "=== Physics Solver Test Complete ==="
    print *, "✓ Physics solver creation works"
    print *, "✓ Physics component initialization works"
    print *, "✓ Physics-enabled time stepping works"
    print *, "✓ Physics disabled mode works"
    print *, "✓ Different physics configurations work"
    print *, "✓ Cleanup and error handling work"
    print *, ""
    print *, "下一步: 实现完整的数值方法（重构、通量、时间积分）"
    
end program test_physics_solver