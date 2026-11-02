! examples/run_eno_weno.f90
program run_eno_weno
    ! Example program: ENO/WENO comparison analysis
    
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction, config_print
    use mesh_module, only: mesh_type
    use registry_module, only: registry_init, registry_cleanup, initialize_default_components
    use physics_solver_module, only: physics_solver, SOLVER_COMPLETED
    
    implicit none
    
    type(cfd_config) :: config_eno3, config_weno3, config_weno5
    type(mesh_type) :: mesh
    type(physics_solver) :: solver_eno3, solver_weno3, solver_weno5
    
    ! Variables to save results
    real(wp) :: time_eno3, time_weno3, time_weno5
    integer :: steps_eno3, steps_weno3, steps_weno5
    integer :: state_eno3, state_weno3, state_weno5
    logical :: all_success
    
    ! Debug: print start marker
    print *, "=========================================="
    print *, "START: ENO/WENO Comparison Analysis"
    print *, "=========================================="
    print *, ""
    
    ! Step 0: Initialize system
    print *, "[STEP 0] Initializing system..."
    print *, "--------------------------------"
    
    call registry_init(verbose=.true.)
    print *, "Registry initialized"
    
    call initialize_default_components()
    print *, "Default components registered"
    print *, ""
    
    ! Step 1: Create mesh
    print *, "[STEP 1] Creating computational mesh..."
    print *, "----------------------------------------"
    
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=40)
    call mesh%print_info()
    print *, ""
    
    ! ========== ENO3 Solver ==========
    print *, "[STEP 2] Configuring and running ENO3 solver..."
    print *, "-----------------------------------------------"
    
    ! Configure ENO3
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
    
    print *, "ENO3 configuration:"
    call config_print(config_eno3)
    print *, ""
    
    ! Create and run ENO3 solver
    print *, "Creating ENO3 solver instance..."
    solver_eno3 = physics_solver(config_eno3, mesh)
    
    print *, "Initializing ENO3 solver..."
    call solver_eno3%initialize()
    
    print *, "Running ENO3 solver..."
    call solver_eno3%run_to_time(config_eno3%final_time)
    
    ! Immediately save ENO3 results
    time_eno3 = solver_eno3%current_time
    steps_eno3 = solver_eno3%current_step
    state_eno3 = solver_eno3%get_state()
    
    print *, "ENO3 solver completed"
    print *, "  Final time: ", time_eno3
    print *, "  Total steps: ", steps_eno3
    print *, "  State: ", state_eno3
    print *, ""
    
    ! ========== WENO3 Solver ==========
    print *, "[STEP 3] Configuring and running WENO3 solver..."
    print *, "------------------------------------------------"
    
    ! Configure WENO3
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
    
    print *, "WENO3 configuration:"
    call config_print(config_weno3)
    print *, ""
    
    ! Create and run WENO3 solver
    print *, "Creating WENO3 solver instance..."
    solver_weno3 = physics_solver(config_weno3, mesh)
    
    print *, "Initializing WENO3 solver..."
    call solver_weno3%initialize()
    
    print *, "Running WENO3 solver..."
    call solver_weno3%run_to_time(config_weno3%final_time)
    
    ! Immediately save WENO3 results
    time_weno3 = solver_weno3%current_time
    steps_weno3 = solver_weno3%current_step
    state_weno3 = solver_weno3%get_state()
    
    print *, "WENO3 solver completed"
    print *, "  Final time: ", time_weno3
    print *, "  Total steps: ", steps_weno3
    print *, "  State: ", state_weno3
    print *, ""
    
    ! ========== WENO5 Solver ==========
    print *, "[STEP 4] Configuring and running WENO5 solver..."
    print *, "------------------------------------------------"
    
    ! Configure WENO5
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
    
    print *, "WENO5 configuration:"
    call config_print(config_weno5)
    print *, ""
    
    ! Create and run WENO5 solver
    print *, "Creating WENO5 solver instance..."
    solver_weno5 = physics_solver(config_weno5, mesh)
    
    print *, "Initializing WENO5 solver..."
    call solver_weno5%initialize()
    
    print *, "Running WENO5 solver..."
    call solver_weno5%run_to_time(config_weno5%final_time)
    
    ! Immediately save WENO5 results
    time_weno5 = solver_weno5%current_time
    steps_weno5 = solver_weno5%current_step
    state_weno5 = solver_weno5%get_state()
    
    print *, "WENO5 solver completed"
    print *, "  Final time: ", time_weno5
    print *, "  Total steps: ", steps_weno5
    print *, "  State: ", state_weno5
    print *, ""
    
    ! ========== Results Summary (BEFORE cleanup!) ==========
    print *, "=========================================="
    print *, "           RESULTS SUMMARY"
    print *, "=========================================="
    print *, ""
    
    print *, "Solver Performance Comparison:"
    print *, "------------------------------"
    
    print *, "ENO3:"
    print *, "  Final time: ", time_eno3
    print *, "  Total steps: ", steps_eno3
    print *, "  State: ", state_eno3
    print *, ""
    
    print *, "WENO3:"
    print *, "  Final time: ", time_weno3
    print *, "  Total steps: ", steps_weno3
    print *, "  State: ", state_weno3
    print *, ""
    
    print *, "WENO5:"
    print *, "  Final time: ", time_weno5
    print *, "  Total steps: ", steps_weno5
    print *, "  State: ", state_weno5
    print *, ""
    
    ! ========== Final Judgment ==========
    print *, "=========================================="
    print *, "           FINAL JUDGMENT"
    print *, "=========================================="
    print *, ""
    
    all_success = (state_eno3 == SOLVER_COMPLETED) .and. &
                  (state_weno3 == SOLVER_COMPLETED) .and. &
                  (state_weno5 == SOLVER_COMPLETED)
    
    if (all_success) then
        print *, "✓ ALL SOLVERS SUCCESSFULLY COMPLETED!"
        print *, ""
        print *, "Parameter Summary:"
        print *, "  Grid cells: ", mesh%ncells
        print *, "  Time step:  ", config_eno3%dt
        print *, "  Final time: ", config_eno3%final_time
        print *, "  Wave speed: ", config_eno3%wave_speed
        print *, "  IC type:    ", trim(config_eno3%ic_type)
        print *, ""
        print *, "Performance Comparison:"
        print *, "  ENO3:  ", steps_eno3, " steps"
        print *, "  WENO3: ", steps_weno3, " steps"
        print *, "  WENO5: ", steps_weno5, " steps"
    else
        print *, "✗ SOME SOLVERS FAILED"
        print *, ""
        print *, "Failure Analysis:"
        if (state_eno3 /= SOLVER_COMPLETED) then
            print *, "  • ENO3 failed with state: ", state_eno3
        end if
        if (state_weno3 /= SOLVER_COMPLETED) then
            print *, "  • WENO3 failed with state: ", state_weno3
        end if
        if (state_weno5 /= SOLVER_COMPLETED) then
            print *, "  • WENO5 failed with state: ", state_weno5
        end if
    end if
    
    print *, ""
    print *, "=========================================="
    print *, "           ANALYSIS COMPLETE"
    print *, "=========================================="
    
    ! ========== Cleanup (AFTER printing results!) ==========
    print *, ""
    print *, "[STEP 5] Cleaning up system..."
    print *, "--------------------------------"
    
    call solver_eno3%cleanup()
    call solver_weno3%cleanup()
    call solver_weno5%cleanup()
    call registry_cleanup()
    
    print *, "All solvers cleaned up"
    print *, "Registry cleaned up"
    print *, ""
    
    ! ========== Wait for user input before exit ==========
    print *, "=========================================="
    print *, "Press ENTER to exit..."
    print *, "=========================================="
    
    ! Wait for user input (uncomment if needed)
    ! read(*,*)
    
    ! Alternative: add a small delay
    print *, "Program will exit in 3 seconds..."
    call sleep(3)
    
    print *, ""
    print *, "=========================================="
    print *, "           PROGRAM END"
    print *, "=========================================="
    
end program run_eno_weno