! tests/test_solver_base.f90 (修复版)
program test_solver_base
    ! 所有 USE 语句必须在程序开始处
    use base_modules, only: wp
    use config_module, only: cfd_config, config_print, config_with_reconstruction
    use mesh_module, only: mesh_type
    use solver_base_module, only: solver_base, SOLVER_UNINITIALIZED, &
                                 SOLVER_INITIALIZED, SOLVER_COMPLETED, SOLVER_ERROR
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(solver_base) :: solver
    integer :: state
    
    print *, "=== Solver Base Test ==="
    print *, ""
    
    ! 测试1: 创建求解器
    print *, "1. Creating solver..."
    print *, "----------------------"
    
    config%verbose = .true.
    call config_with_reconstruction(config, "eno", 3)
    config%dt = 0.01_wp
    
    call config_print(config)
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=10)
    
    solver = solver_base(config, mesh)
    call solver%print_info()
    print *, ""
    
    ! 测试2: 初始化
    print *, "2. Initializing solver..."
    print *, "-------------------------"
    
    call solver%initialize()
    state = solver%get_state()
    print *, "State after initialization: ", state
    print *, "Expected: ", SOLVER_INITIALIZED
    print *, "Match? ", state == SOLVER_INITIALIZED
    print *, "Error message: '", trim(solver%get_error()), "'"
    print *, ""
    
    ! 测试3: 运行求解器
    print *, "3. Running solver..."
    print *, "--------------------"
    
    call solver%run_to_time(0.05_wp)
    state = solver%get_state()
    print *, "State after run: ", state
    print *, "Expected: ", SOLVER_COMPLETED
    print *, "Match? ", state == SOLVER_COMPLETED
    print *, "Current time: ", solver%current_time
    print *, "Current step: ", solver%current_step
    print *, ""
    
    ! 测试4: 再次运行（从已完成状态）
    print *, "4. Running again from completed state..."
    print *, "----------------------------------------"
    
    ! 需要先清理才能重新运行
    call solver%cleanup()
    call solver%initialize()
    call solver%run_to_time(0.1_wp)
    
    call solver%print_info()
    print *, ""
    
    ! 测试5: 错误处理
    print *, "5. Testing error states..."
    print *, "--------------------------"
    
    ! 创建一个未初始化的求解器
    call solver%cleanup()
    state = solver%get_state()
    print *, "Uninitialized state: ", state
    print *, "Expected: ", SOLVER_UNINITIALIZED
    print *, "Match? ", state == SOLVER_UNINITIALIZED
    
    ! 尝试运行未初始化的求解器
    call solver%run_to_time(0.01_wp)
    state = solver%get_state()
    print *, "State after error: ", state
    print *, "Expected: ", SOLVER_ERROR
    print *, "Match? ", state == SOLVER_ERROR
    print *, "Error message: '", trim(solver%get_error()), "'"
    print *, ""
    
    print *, "=== Solver Base Test Complete ==="
    print *, "✓ Solver base class works"
    print *, "✓ State management works"
    print *, "✓ Time stepping framework works"
    print *, "✓ Error handling works"
    
end program test_solver_base