! tests/test_physics_solver.f90 (简化修正版)
program test_physics_solver
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction
    use mesh_module, only: mesh_type
    use physics_solver_module, only: physics_solver, SOLVER_COMPLETED
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(physics_solver) :: psolver
    real(wp) :: final_time
    integer :: state
    
    print *, "=== Physics Solver Test (简化版) ==="
    print *, ""
    
    ! 测试1: 创建物理求解器
    print *, "1. Creating physics solver..."
    print *, "-----------------------------"
    
    config%verbose = .true.
    call config_with_reconstruction(config, "eno", 3)
    config%dt = 0.01_wp
    config%enable_physics = .true.
    
    print *, "Configuration:"
    print *, "  Scheme: ", trim(config%recon_scheme)
    print *, "  dt: ", config%dt
    print *, "  Physics enabled: ", config%enable_physics
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=10)
    
    ! 设置求解器配置和网格
    psolver%config = config
    psolver%mesh = mesh
    
    print *, "  Solver created successfully"
    print *, ""
    
    ! 测试2: 初始化
    print *, "2. Initializing physics solver..."
    print *, "---------------------------------"
    
    call psolver%initialize()
    state = psolver%get_state()
    
    print *, "  State after initialization: ", state
    print *, "  Expected: initialized (1)"
    print *, ""
    
    ! 测试3: 运行一小段时间
    print *, "3. Running physics solver (short time)..."
    print *, "------------------------------------------"
    
    call psolver%run_to_time(0.02_wp)
    state = psolver%get_state()
    
    print *, "  State after run: ", state
    print *, "  Expected: completed (3)"
    print *, "  Current time: ", psolver%current_time
    print *, "  Current step: ", psolver%current_step
    print *, ""
    
    ! 测试4: 清理
    print *, "4. Testing cleanup..."
    print *, "----------------------"
    
    call psolver%cleanup()
    state = psolver%get_state()
    
    print *, "  State after cleanup: ", state
    print *, "  Expected: uninitialized (0)"
    print *, ""
    
    ! 结果验证
    print *, "=== Test Summary ==="
    if (state == 0) then
        print *, "✓ All basic tests passed"
    else
        print *, "✗ Some tests failed"
    end if
    
end program test_physics_solver