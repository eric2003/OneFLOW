! ==================== cfd_solver.f90 ====================
! OneFLOW-CFD Solver for 1D Convection Equation

module cfd_solver
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType, update_oldfield
    use boundary_module, only: boundary
    use initial_condition_module, only: initial_condition, analytical_solution, init_ic
	use recon_module, only: ReconType
    use time_integration_module, only: runge_kutta
    use string_utils, only: to_lower, to_upper
	use io_utils, only: write_simulation_results
    
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
        
        call init_ic(config%ic_type)
        
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
        integer :: ist, ied
        integer :: n
        real(dp), allocatable :: um(:)
        
        ist = this%domain%ist
        ied = this%domain%ied
        n = ied - ist
        allocate(um(n))
        
        do i = ist, ied - 1
            j = i - ist + 1
            this%solution%u(i) = initial_condition(this%domain%mesh%xcc(j))
            um(j) = this%solution%u(i)
        end do
        
       
        call save_results(this%domain%mesh%xcc, um, n, "precise_wave.dat")
        

        ! 使用 boundary_module 中的函数
        call boundary(this%solution%u, this%domain%nghosts, this%domain%ist, this%domain%ied)
        call update_oldfield(this%solution%un, this%solution%u, this%domain%ntcells)
    end subroutine init_field
    
   subroutine save_results(x, u, n, filename)
        real(dp), intent(in) :: x(:), u(:)
        integer, intent(in) :: n
        character(len=*), intent(in) :: filename
        integer :: i
        
        open(20, file=filename, status='replace')
        write(20, *) "# Complex wave with precise parameters"
        write(20, *) "# a=0.5, z=-0.7, δ=0.005, α=10, β=log2/(36δ²)"
        write(20, *) "# x u"
        do i = 1, n
            write(20, '(2ES16.8)') x(i), u(i)
        end do
        close(20)
        print *, "  Saved to ", trim(filename)
    end subroutine        
    
     
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
        character(len=:), allocatable :: uscheme
		character(len=:), allocatable :: lscheme
		character(len=:), allocatable :: filename  ! 自动管理长度

        allocate(u_result(this%domain%mesh%ncells))
        
        call this%init_field()
        
        config => this%config
        
        uscheme = to_upper(this%config%recon_scheme)
        
        ! Run simulation
        print *, "=========================================="
        print *, "Running " // trim(uscheme) // " simulation..."
        print *, "  Scheme:" // trim(uscheme), config%spatial_order
        print *, "  Time step: ", config%dt
        print *, "=========================================="        
        
        final_time = config%final_time
        
        t = 0.0_dp        
        dt_old = config%dt
        dt = dt_old
        max_steps = 10000000  ! Safety limit
        
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
		lscheme = to_upper(this%config%recon_scheme)
		filename = trim(lscheme) // "_results.txt"
		
		print *, "Writing results to files..."
		
		call write_simulation_results(filename, this%domain%mesh%xcc, u_result, trim(uscheme))
		
    end function run_simulation
    
end module cfd_solver