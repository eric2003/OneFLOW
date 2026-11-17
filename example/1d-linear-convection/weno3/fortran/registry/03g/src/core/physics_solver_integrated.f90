! src/solver/physics_solver_integrated.f90
module physics_solver_integrated
    use base_modules, only: wp, ip
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create
    use registry_module, only: create_component
    
    implicit none
    private
    
    type, public :: physics_solver_integrated
        ! 核心组件
        type(cfd_config) :: config
        type(mesh_type) :: mesh
        type(domain_type) :: domain
        type(solution_type) :: solution
        
        ! 物理组件
        class(*), allocatable :: equation
        class(*), allocatable :: problem
        
        ! 数值组件
        class(*), allocatable :: reconstructor
        class(*), allocatable :: flux_calculator
        class(*), allocatable :: boundary_condition
        class(*), allocatable :: residual_calculator
        
    contains
        procedure :: initialize => solver_initialize
        procedure :: run => solver_run
        procedure :: cleanup => solver_cleanup
        procedure, private :: create_components
        procedure, private :: apply_boundary
    end type
    
contains

    subroutine solver_initialize(this)
        class(physics_solver_integrated), intent(inout) :: this
        
        ! 创建域和解
        this%domain = domain_create(this%config, this%mesh)
        this%solution = solution_create(this%domain)
        
        ! 创建所有组件
        call this%create_components()
        
        ! 应用初始条件
        call this%apply_initial_condition()
        
        print *, "[SOLVER] Initialized with physics integration"
    end subroutine
    
    subroutine create_components(this)
        class(physics_solver_integrated), intent(inout) :: this
        
        ! 创建方程
        call create_component("equation", "linear_advection", this%config, this%equation)
        
        ! 创建问题
        call create_component("problem", "linear_advection", this%config, this%problem)
        
        ! 创建数值组件
        call create_component("reconstructor", this%config%recon_scheme, &
                              this%config, this%reconstructor)
        call create_component("flux", this%config%flux_type, &
                              this%config, this%flux_calculator)
        call create_component("boundary", this%config%boundary_type, &
                              this%cfd_context(), this%boundary_condition)
    end subroutine
    
    subroutine apply_initial_condition(this)
        class(physics_solver_integrated), intent(inout) :: this
        
        ! 通过问题创建初始条件
        select type(prob => this%problem)
        type is (linear_advection_problem)
            class(*), allocatable :: ic
            call prob%create_ic(this%config, ic)
            
            ! 应用初始条件到解
            select type(ic_inst => ic)
            type is (step_function_ic)
                call ic_inst%apply(this%solution)
            end select
        end select
    end subroutine
    
    function cfd_context(this) result(ctx)
        class(physics_solver_integrated), intent(in) :: this
        type(cfd_context_type) :: ctx
        
        ! 创建包含求解器所有组件的上下文
        ctx%config => this%config
        ctx%domain => this%domain
        ctx%solution => this%solution
        ctx%equation => this%equation
    end function
    
    subroutine solver_run(this, final_time)
        class(physics_solver_integrated), intent(inout) :: this
        real(wp), intent(in) :: final_time
        
        real(wp) :: t, dt, dt_original
        integer :: step
        
        dt_original = this%config%dt
        t = 0.0_wp
        step = 0
        
        do while (t < final_time - 1e-12_wp)
            dt = min(this%config%dt, final_time - t)
            
            ! 时间步进（需要实现）
            ! call this%time_step(dt)
            
            t = t + dt
            step = step + 1
        end do
        
        this%config%dt = dt_original
        print *, "[SOLVER] Completed at t = ", t, ", steps = ", step
    end subroutine
    
end module