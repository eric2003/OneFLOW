! Perform ENO-WENO comparative analysis
module config_factory_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    implicit none
    
    private
    
    public :: create_eno_config, create_weno_config
    
contains
    
    subroutine create_eno_config(config)
        type(CfdConfigType), intent(out) :: config
        
        config%recon_scheme = "eno"
        config%spatial_order = 3
        config%flux_type = 0
        config%rk_order = 2
        config%wave_speed = 1.0_dp
        config%final_time = 0.625_dp
        config%cfl = 1.0_dp
        config%dt = 0.0025_dp
    end subroutine create_eno_config
    
    subroutine create_weno_config(config)
        type(CfdConfigType), intent(out) :: config
        
        config%recon_scheme = "weno"
        config%spatial_order = 3
        config%flux_type = 0
        config%rk_order = 2
        config%wave_speed = 1.0_dp
        config%final_time = 0.625_dp
        config%cfl = 1.0_dp
        config%dt = 0.0025_dp
    end subroutine create_weno_config
    
end module config_factory_module
    
! Perform ENO-WENO comparative analysis
subroutine performEnoWenoAnalysis()
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    use domain_module, only: ComputationalDomainType
    use initial_condition_module, only: analytical_solution
    use cfd_solver, only: CfdType, init_field, run_simulation
    use config_factory_module, only: create_eno_config, create_weno_config
    
    type(CfdConfigType) :: config_eno3, config_weno3
    type(MeshType), target :: mesh
    type(CfdType) :: cfd_eno3, cfd_weno3
    real(dp), allocatable :: u_eno(:), u_weno(:), u_analytical(:)
    real(dp), pointer :: xcc(:)
    integer :: i, ncells, iunit
        
    ! Initialize mesh
    call mesh%init()
    ncells = mesh%ncells
    xcc => mesh%xcc
    
        
    print *, "=========================================="
    print *, "Mesh parameters:"
    print *, "  ncells = ", ncells
    print *, "  dx = ", mesh%dx
    print *, "  L = ", mesh%L
    print *, "=========================================="
        
    ! Configure ENO3
	call create_eno_config(config_eno3)
    
    ! Configure WENO3
    call create_weno_config(config_weno3)
        
    ! Create CFD solvers
    call cfd_eno3%init(config_eno3, mesh)
    call cfd_weno3%init(config_weno3, mesh)
        
    ! Allocate arrays
    allocate(u_eno(ncells), u_weno(ncells), u_analytical(ncells))
        
    ! Run ENO simulation
    print *, "=========================================="
    print *, "Running ENO simulation..."
    print *, "  Scheme: ENO", config_eno3%spatial_order
    print *, "  Time step: ", config_eno3%dt
    print *, "=========================================="
        
    !call init_field(cfd_eno3)
    call cfd_eno3.my_init_field()
    
    u_eno = run_simulation(cfd_eno3, config_eno3%final_time)
        
    ! Run WENO simulation
    print *, "=========================================="
    print *, "Running WENO simulation..."
    print *, "  Scheme: WENO", config_weno3%spatial_order
    print *, "  Time step: ", config_weno3%dt
    print *, "=========================================="
        
    call init_field(cfd_weno3)
    u_weno = run_simulation(cfd_weno3, config_weno3%final_time)
        
    ! Compute analytical solution
    print *, "Computing analytical solution..."
    do i = 1, ncells
        u_analytical(i) = analytical_solution(xcc(i), config_weno3%final_time, &
                                            config_weno3%wave_speed, mesh%L)
    end do
        
    ! Write results to files
    print *, "Writing results to files..."
        
    ! Write ENO results
    open(newunit=iunit, file='eno_results.txt', status='replace')
    write(iunit, '(A)') '# x, u (ENO)'
    do i = 1, ncells
        write(iunit, '(2F12.6)') xcc(i), u_eno(i)
    end do
    close(iunit)
        
    ! Write WENO results
    open(newunit=iunit, file='weno_results.txt', status='replace')
    write(iunit, '(A)') '# x, u (WENO)'
    do i = 1, ncells
        write(iunit, '(2F12.6)') xcc(i), u_weno(i)
    end do
    close(iunit)
        
    ! Write analytical results
    open(newunit=iunit, file='analytical_results.txt', status='replace')
    write(iunit, '(A)') '# x, u (Analytical)'
    do i = 1, ncells
        write(iunit, '(2F12.6)') xcc(i), u_analytical(i)
    end do
    close(iunit)
        
    ! Print some statistics
    print *, "=========================================="
    print *, "Simulation statistics:"
    print *, "  ENO min/max: ", minval(u_eno), maxval(u_eno)
    print *, "  WENO min/max: ", minval(u_weno), maxval(u_weno)
    print *, "  Analytical min/max: ", minval(u_analytical), maxval(u_analytical)
    print *, "=========================================="
        
    print *, "Simulation completed successfully!"
    print *, "Results written to:"
    print *, "  eno_results.txt"
    print *, "  weno_results.txt"
    print *, "  analytical_results.txt"
    print *, ""
    print *, "To generate the comparison plot, run:"
    print *, "  python postprocess.py"
    print *, "=========================================="
        
    deallocate(u_eno, u_weno, u_analytical)
end subroutine performEnoWenoAnalysis
    
! ===================================================================
! Main Program
! ===================================================================
program main
    implicit none
    
    print *, "=========================================="
    print *, "OneFLOW-CFD Solver for 1D Convection"
    print *, "ENO vs WENO Comparison"
    print *, "=========================================="
    
    call performEnoWenoAnalysis()
    
    print *, "Program finished successfully!"
    
end program main