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
        type(CfdConfigType), pointer :: config => null() 
        type(ComputationalDomainType) :: domain
        type(SolutionType) :: solution
        type(MeshType)  :: mesh
        type(ReconType), allocatable :: recon
    contains
        procedure :: init => cfd_init
        procedure :: init_field => init_field
        procedure :: run => run_simulation
    end type CfdType
    
    ! ===================================================================
    ! Public Procedures
    ! ===================================================================
    
    contains
    
    ! ===================================================================
    ! Initialization Methods
    ! ===================================================================
    
    ! CFD solver initialization
   
    subroutine cfd_init(this, config, mesh)
        class(CfdType), intent(inout) :: this
        type(CfdConfigType), intent(in), target :: config
        type(MeshType), intent(in) :: mesh
        
        this%config => config
        this%mesh = mesh
        call this%domain%init(config, mesh)
        
        call this%solution%init(this%domain)
		
        allocate(this%recon)
        call this%recon%init(config, this%domain)
        
        ! Adjust time step based on CFL condition
        call calculate_dt(this)
    end subroutine cfd_init    
    
    ! Calculate time step based on CFL condition
    subroutine calculate_dt(cfd)
        type(CfdType), intent(inout) :: cfd
        type(CfdConfigType), pointer :: config
        real(dp) :: dt_cfl
        
        config => cfd%config
        
        ! CFL condition: dt <= CFL * dx / |wave_speed|
        dt_cfl = config%cfl * cfd%domain%mesh%dx / abs(config%wave_speed)
        
        if (config%dt > dt_cfl) then
            print *, "Adjusting time step for stability:"
            print *, "  Original dt = ", config%dt
            print *, "  CFL dt = ", dt_cfl
            config%dt = dt_cfl
            print *, "  Using dt = ", config%dt
        end if

    end subroutine calculate_dt
    
    ! ===================================================================
    ! Field Initialization (使用外部模块的函数)
    ! ===================================================================
    
    ! Initialize field with step function
    subroutine init_field(this)
        class(CfdType), intent(inout) :: this
        integer :: i, j
        integer ::ist, ied
        ist = this%domain%ist
        ied = this%domain%ied
        
        do i = ist, ied - 1
            j = i - ist + 1
            this%solution%u(i) = initial_condition(this%domain%mesh%xcc(j))
        end do
        

        ! 使用 boundary_module 中的函数
        call boundary(this%solution%u, this%domain%nghosts, this%domain%ist, this%domain%ied)
        call update_oldfield(this%solution%un, this%solution%u, this%domain%ntcells)
    end subroutine init_field
    
     
    ! ===================================================================
    ! Simulation Driver
    ! ===================================================================
    
    ! Run simulation to final time
    function run_simulation(this) result(u_result)
        class(CfdType), intent(inout) :: this
        
        real(dp), allocatable :: u_result(:)
        
        real(dp) :: t, dt, dt_old, final_time
        integer :: step, max_steps
        
        type(CfdConfigType), pointer :: config
        
        
        allocate(u_result(this%domain%mesh%ncells))
        
        config => this%config
        final_time = config%final_time
        
        t = 0.0_dp        
        dt_old = config%dt
        dt = dt_old
        max_steps = 10000  ! Safety limit
        
        print *, "Starting time integration..."
        print *, "  Final time: ", final_time
        print *, "  Time step: ", dt
        print *, "  CFL number: ", config%cfl
        
        step = 0
        do while (t < final_time - 1.0e-12_dp .and. step < max_steps)
            step = step + 1
            
            if (t + dt > final_time) then
                dt = final_time - t
            end if
            
            config%dt = dt
            call runge_kutta(this%recon, config, this%domain, this%solution)
            t = t + dt
            
            ! Progress report
            if (mod(step, 100) == 0) then
                print *, "  Step ", step, ", Time = ", t
            end if
        end do
        
        if (step >= max_steps) then
            print *, "Warning: Reached maximum number of steps (", max_steps, ")"
        end if
        
        config%dt = dt_old
        
        print *, "Time integration completed:"
        print *, "  Total steps: ", step
        print *, "  Final time: ", t
        
        ! Extract physical solution (without ghost cells)
        u_result = this%solution%u(this%domain%ist:this%domain%ied)
    end function run_simulation
    
end module cfd_solver