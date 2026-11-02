! ==================== cfd_solver.f90 ====================
! OneFLOW-CFD Solver for 1D Convection Equation

module cfd_solver
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType, update_oldfield
    use boundary_module, only: boundary
    use initial_condition_module, only: initial_condition, analytical_solution
	use recon_module, only: init_coef, ReconType,  reconstruction
    use time_integration_module, only: runge_kutta
    
    implicit none
    
    private

    ! ===================================================================
    ! Main CFD Solver Class
    ! ===================================================================
    type, public :: CfdType
        !type(CfdConfigType) :: config
        type(CfdConfigType), pointer :: config => null() 
        type(ComputationalDomainType) :: domain
        type(SolutionType) :: solution
        type(ReconType), allocatable :: recon
    contains
        procedure :: init => cfd_init
    end type CfdType
    
    ! ===================================================================
    ! Public Procedures
    ! ===================================================================
    public :: run_simulation, init_field
    
    contains
    
    ! ===================================================================
    ! Initialization Methods
    ! ===================================================================
    
    ! CFD solver initialization
   
    subroutine cfd_init(this, config, domain)
        class(CfdType), intent(inout) :: this
        type(CfdConfigType), intent(in), target :: config
        type(ComputationalDomainType), intent(in) :: domain
        
        this%config => config
        this%domain = domain
        call this%solution%init(domain)
		
        allocate(this%recon)
        call this%recon%init(config, domain)
        
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
        integer ::ist, ied
        ist = cfd%domain%ist
        ied = cfd%domain%ied
        
        do i = ist, ied - 1
            j = i - ist + 1
            cfd%solution%u(i) = initial_condition(cfd%domain%mesh%xcc(j))
        end do
        

        ! 使用 boundary_module 中的函数
        call boundary(cfd%solution%u, cfd%domain%nghosts, cfd%domain%ist, cfd%domain%ied)
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine init_field
    
     
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
        
        type(CfdConfigType), pointer :: config
        
        
        allocate(u_result(cfd%domain%mesh%ncells))
        
        t = 0.0_dp
        
        !print *, "cfd%config : ", loc(cfd%config)

        config => cfd%config
        !dt_old = cfd%config%dt
        dt_old = config%dt
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
            call runge_kutta(cfd%recon, cfd%config, cfd%domain, cfd%solution)
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
    
end module cfd_solver