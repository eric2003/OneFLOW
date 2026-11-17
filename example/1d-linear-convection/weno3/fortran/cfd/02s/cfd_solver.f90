! ==================== cfd_solver.f90 ====================
! OneFLOW-CFD Solver for 1D Convection Equation

module cfd_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use config_module, only: CfdConfigType, create_eno_config, create_weno_config
    use mesh_module, only: MeshType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    use flux_module, only: inviscid_flux
    use boundary_module, only: boundary
    use initial_condition_module, only: initial_condition, analytical_solution
	use recon_module, only: init_coef, ReconType,  &
        eno_reconstruct, weno_reconstruct
    implicit none
    
    private

    ! ===================================================================
    ! Main CFD Solver Class
    ! ===================================================================
    ! ======================================================
    ! 【修改3】：CfdType 中去掉多态 reconstructor，改为两个独立变量
    ! ======================================================
    type, public :: CfdType
        type(CfdConfigType) :: config
        type(ComputationalDomainType) :: domain
        type(SolutionType) :: solution
        type(ReconType), allocatable :: my_recon
    contains
        procedure :: init => cfd_init
    end type CfdType
	
    
    ! ===================================================================
    ! Public Procedures
    ! ===================================================================
    public :: run_simulation, init_field, performEnoWenoAnalysis
    
    
    contains
    
    ! ===================================================================
    ! Initialization Methods
    ! ===================================================================
    
    ! CFD solver initialization
   
    subroutine cfd_init(this, config, domain)
        class(CfdType), intent(inout) :: this
        type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        
        this%config = config
        this%domain = domain
        call this%solution%init(domain)
		
        allocate(this%my_recon)
        call this%my_recon%init(config, domain)
        
        ! Adjust time step based on CFL condition
        call calculate_dt(this)
    end subroutine cfd_init    
    
    ! Calculate time step based on CFL condition
    subroutine calculate_dt(cfd)
        type(CfdType), intent(inout) :: cfd
        
        real(dp) :: dt_cfl
        
        ! CFL condition: dt <= CFL * dx / |wave_speed|
        dt_cfl = cfd%config%cfl * cfd%domain%mesh%dx / abs(cfd%config%wave_speed)
        
        if (cfd%config%dt > dt_cfl) then
            print *, "Adjusting time step for stability:"
            print *, "  Original dt = ", cfd%config%dt
            print *, "  CFL dt = ", dt_cfl
            cfd%config%dt = dt_cfl
            print *, "  Using dt = ", cfd%config%dt
        end if

    end subroutine calculate_dt
    
    ! ===================================================================
    ! Field Initialization (使用外部模块的函数)
    ! ===================================================================
    
    ! Initialize field with step function
    subroutine init_field(cfd)
        type(CfdType), intent(inout) :: cfd
        integer :: i, j
        
        do i = 1, cfd%domain%mesh%ncells
            ! 使用 initial_condition_module 中的函数
            cfd%solution%u(cfd%domain%ist + i - 1) = initial_condition(cfd%domain%mesh%xcc(i))
        end do
        
        ! 使用 boundary_module 中的函数
        call boundary(cfd%solution%u, cfd%domain%nghosts, cfd%domain%ist, cfd%domain%ied)
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine init_field
    
    ! ===================================================================
    ! Reconstruction Methods
    ! ===================================================================
    
    subroutine reconstruction(q, recon, config, domain, solution)
        real(dp), intent(in) :: q(:)
		type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        type(ReconType), intent(inout) :: recon
		
		call recon%reconstruct(q, config, domain, solution)
    end subroutine reconstruction
    
    ! ===================================================================
    ! Residual Computation (更新：使用新的通量模块)
    ! ===================================================================
    
    ! Compute residual (flux divergence)
    subroutine residual(q, cfd)
        real(dp), intent(in) :: q(:)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i
        
        ! Reconstruction using current solution
        call reconstruction(q, cfd%my_recon, cfd%config, cfd%domain, cfd%solution)
        
        ! Compute fluxes - 使用新的通量模块
        call inviscid_flux(cfd%solution%q_face_left, cfd%solution%q_face_right, &
                          cfd%solution%flux, cfd%config, cfd%domain%mesh)
        
        ! Compute residual
        do i = 1, cfd%domain%mesh%ncells
            cfd%solution%res(i) = -(cfd%solution%flux(i+1) - cfd%solution%flux(i)) / &
                                   cfd%domain%mesh%dx
        end do
    end subroutine residual
    
    ! ===================================================================
    ! Time Integration (保持原有实现)
    ! ===================================================================
    
    ! Update old field
    subroutine update_oldfield(qn, q, n)
        real(dp), intent(out) :: qn(:)
        real(dp), intent(in) :: q(:)
        integer, intent(in) :: n
        
        qn(1:n) = q(1:n)
    end subroutine update_oldfield
    
    ! 1st-order Runge-Kutta (Euler)
    subroutine runge_kutta_1(cfd)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i, j
        real(dp) :: dt
        
        dt = cfd%config%dt
        
        ! Compute residual
        call residual(cfd%solution%u, cfd)
        
        ! Update solution
        do i = cfd%domain%ist, cfd%domain%ied - 1
            j = i - cfd%domain%ist + 1
            cfd%solution%u(i) = cfd%solution%u(i) + dt * cfd%solution%res(j)
        end do
        
        ! Apply boundary conditions again
        call boundary(cfd%solution%u, cfd%domain%nghosts, cfd%domain%ist, cfd%domain%ied)
        
        ! Save old solution
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine runge_kutta_1
    
    ! 2nd-order Runge-Kutta (Heun)
    subroutine runge_kutta_2(cfd)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i, j
        real(dp) :: dt
        
        dt = cfd%config%dt
        
        ! Stage 1
        call residual(cfd%solution%u, cfd)
        
        do i = cfd%domain%ist, cfd%domain%ied - 1
            j = i - cfd%domain%ist + 1
            cfd%solution%u(i) = cfd%solution%u(i) + dt * cfd%solution%res(j)
        end do
        call boundary(cfd%solution%u, cfd%domain%nghosts, cfd%domain%ist, cfd%domain%ied)
        
        ! Stage 2
        call residual(cfd%solution%u, cfd)
        do i = cfd%domain%ist, cfd%domain%ied - 1
            j = i - cfd%domain%ist + 1
            cfd%solution%u(i) = 0.5_dp * cfd%solution%un(i) + &
                                0.5_dp * cfd%solution%u(i) + &
                                0.5_dp * dt * cfd%solution%res(j)
        end do
        call boundary(cfd%solution%u, cfd%domain%nghosts, cfd%domain%ist, cfd%domain%ied)
        
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine runge_kutta_2
    
    ! Runge-Kutta selection
    subroutine runge_kutta(cfd)
        type(CfdType), intent(inout) :: cfd
        
        select case(cfd%config%rk_order)
        case(1)
            call runge_kutta_1(cfd)
        case(2)
            call runge_kutta_2(cfd)
        case default
            call runge_kutta_1(cfd)
        end select
    end subroutine runge_kutta
    
    ! ===================================================================
    ! Simulation Driver
    ! ===================================================================
    
    ! Run simulation to final time
    function run_simulation(cfd, final_time) result(u_result)
        type(CfdType), intent(inout) :: cfd
        real(dp), intent(in) :: final_time
        real(dp), allocatable :: u_result(:)
        
        real(dp) :: t, dt, dt_old
        integer :: step, max_steps
        
        allocate(u_result(cfd%domain%mesh%ncells))
        
        t = 0.0_dp
        dt_old = cfd%config%dt
        dt = dt_old
        max_steps = 10000  ! Safety limit
        
        print *, "Starting time integration..."
        print *, "  Final time: ", final_time
        print *, "  Time step: ", dt
        print *, "  CFL number: ", cfd%config%cfl
        
        step = 0
        do while (t < final_time - 1.0e-12_dp .and. step < max_steps)
            step = step + 1
            
            if (t + dt > final_time) then
                dt = final_time - t
            end if
            
            cfd%config%dt = dt
            call runge_kutta(cfd)
            t = t + dt
            
            ! Progress report
            if (mod(step, 100) == 0) then
                print *, "  Step ", step, ", Time = ", t
            end if
        end do
        
        if (step >= max_steps) then
            print *, "Warning: Reached maximum number of steps (", max_steps, ")"
        end if
        
        cfd%config%dt = dt_old
        
        print *, "Time integration completed:"
        print *, "  Total steps: ", step
        print *, "  Final time: ", t
        
        ! Extract physical solution (without ghost cells)
        u_result = cfd%solution%u(cfd%domain%ist:cfd%domain%ied)
    end function run_simulation
    
    ! ===================================================================
    ! Main Analysis Function
    ! ===================================================================
    
    ! Perform ENO-WENO comparative analysis
    subroutine performEnoWenoAnalysis()
        type(CfdConfigType) :: config_eno3, config_weno3
        type(MeshType) :: mesh
        type(ComputationalDomainType) :: domain_eno3, domain_weno3
        type(CfdType) :: cfd_eno3, cfd_weno3
        real(dp), allocatable :: u_eno(:), u_weno(:), u_analytical(:)
        real(dp), allocatable :: xcc(:)
        integer :: i, ncells, iunit
        
        ! Initialize mesh
        call mesh%init()
        ncells = mesh%ncells
        allocate(xcc(ncells))
        xcc = mesh%xcc
        
        print *, "=========================================="
        print *, "Mesh parameters:"
        print *, "  ncells = ", ncells
        print *, "  dx = ", mesh%dx
        print *, "  L = ", mesh%L
        print *, "=========================================="
        
        ! Configure ENO3
        config_eno3 = create_eno_config()
        
        ! Configure WENO3
        config_weno3 = create_weno_config()
        
        ! Create domains
        call domain_eno3%init(mesh, config_eno3)
        call domain_weno3%init(mesh, config_weno3)
        
        ! Create CFD solvers
        call cfd_eno3%init(config_eno3, domain_eno3)
        call cfd_weno3%init(config_weno3, domain_weno3)
        
        ! Allocate arrays
        allocate(u_eno(ncells), u_weno(ncells), u_analytical(ncells))
        
        ! Run ENO simulation
        print *, "=========================================="
        print *, "Running ENO simulation..."
        print *, "  Scheme: ENO", config_eno3%spatial_order
        print *, "  Time step: ", config_eno3%dt
        print *, "=========================================="
        
        call init_field(cfd_eno3)
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
        
        deallocate(u_eno, u_weno, u_analytical, xcc)
    end subroutine performEnoWenoAnalysis

end module cfd_solver

! ===================================================================
! Main Program
! ===================================================================
program main
    use cfd_solver
    implicit none
    
    print *, "=========================================="
    print *, "OneFLOW-CFD Solver for 1D Convection"
    print *, "ENO vs WENO Comparison"
    print *, "=========================================="
    
    call performEnoWenoAnalysis()
    
    print *, "Program finished successfully!"
    
end program main