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
        !config%final_time = 0.625_dp
        config%final_time = 8.0_dp
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
        !config%final_time = 0.625_dp
        config%final_time = 8.0_dp
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
    use cfd_solver, only: CfdType
    use config_factory_module, only: create_eno_config, create_weno_config
	use io_utils, only: write_simulation_results
    
    type(CfdConfigType) :: config_eno3, config_weno3
    type(MeshType), target :: mesh
    type(CfdType) :: cfd_eno3, cfd_weno3
    real(dp), allocatable :: u_eno(:), u_weno(:), u_analytical(:)
    real(dp), pointer :: xcc(:)
    integer :: i, ncells, iunit
        
    ! Initialize mesh
    !call mesh%init()
    call mesh%init(-1.0_dp, 1.0_dp, 200)
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
    call cfd_eno3%init(config_eno3, mesh)
    
    ! Configure WENO3
    call create_weno_config(config_weno3)
    call cfd_weno3%init(config_weno3, mesh)
        
    ! Allocate arrays
    allocate(u_eno(ncells), u_weno(ncells), u_analytical(ncells))
        
    u_eno = cfd_eno3.run()
         
    u_weno = cfd_weno3.run()
        
    ! Compute analytical solution
    print *, "Computing analytical solution..."
    do i = 1, ncells
        u_analytical(i) = analytical_solution(xcc(i), config_weno3%final_time, &
                                            config_weno3%wave_speed, mesh%xmin, mesh%xmax)
    end do
        
    ! Write results to files
    print *, "Writing results to files..."
	
	call write_simulation_results('analytical_results.txt', mesh%xcc, u_analytical, 'Analytical')
        
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