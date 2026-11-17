! src/solver/physics_solver.f90 (简化版)
module physics_solver_module
    use base_modules, only: wp => wp, ip => ip
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create
    
    use physics_interface, only: physics_equation, physics_problem
    use linear_convection_equation, only: linear_convection_eq, create_linear_convection_eq
    use linear_convection_problem, only: linear_convection_prob, create_linear_convection_prob
    
    use component_manager_module, only: create_reconstructor, create_flux_calculator
    use reconstructor_base_module, only: reconstructor_base
    use flux_base_module, only: flux_calculator_base
    
    implicit none
    private
    
    ! 明确导出列表
    public :: wp, ip, physics_solver, create_physics_solver
    public :: SOLVER_UNINITIALIZED, SOLVER_INITIALIZED, SOLVER_RUNNING
    public :: SOLVER_COMPLETED, SOLVER_ERROR
    
    ! 求解器状态枚举
    integer, parameter :: SOLVER_UNINITIALIZED = 0
    integer, parameter :: SOLVER_INITIALIZED = 1
    integer, parameter :: SOLVER_RUNNING = 2
    integer, parameter :: SOLVER_COMPLETED = 3
    integer, parameter :: SOLVER_ERROR = 4
    
    ! 物理求解器类型（不继承，独立实现）
    type :: physics_solver
        ! 基本组件
        type(cfd_config) :: config
        type(mesh_type) :: mesh
        type(domain_type) :: domain
        type(solution_type) :: solution
        
        ! 物理组件
        class(physics_equation), allocatable :: equation
        class(physics_problem), allocatable :: problem
        
        ! 数值组件
        class(reconstructor_base), allocatable :: reconstructor
        class(flux_calculator_base), allocatable :: flux_calculator
        
        ! 状态管理
        integer :: state = SOLVER_UNINITIALIZED
        character(len=100) :: error_message = ""
        real(wp) :: current_time = 0.0_wp
        integer(ip) :: current_step = 0
        
        ! 时间控制
        real(wp) :: dt_original = 0.0_wp
        logical :: physics_initialized = .false.
    contains
        procedure :: initialize => physics_solver_initialize
        procedure :: step => physics_solver_step
        procedure :: run_to_time => physics_solver_run_to_time
        procedure :: cleanup => physics_solver_cleanup
        procedure :: get_state => physics_solver_get_state
        procedure :: get_error => physics_solver_get_error
        procedure :: print_info => physics_solver_print_info
        procedure, private :: create_physics_components
        procedure, private :: create_numerical_components
    end type physics_solver
    
    ! 构造函数接口
    interface physics_solver
        module procedure create_physics_solver
    end interface
    
contains

    ! ==================== 构造函数 ====================
    
    function create_physics_solver(config, mesh) result(solver)
        type(cfd_config), intent(in) :: config
        type(mesh_type), intent(in) :: mesh
        type(physics_solver) :: solver
        
        solver%config = config
        solver%mesh = mesh
        
        ! 创建域
        solver%domain = domain_create(config, mesh)
        
        ! 创建解
        solver%solution = solution_create(solver%domain)
        
        ! 保存原始时间步长
        solver%dt_original = config%dt
        
        ! 创建组件
        call solver%create_physics_components()
        call solver%create_numerical_components()
        
        if (config%verbose) then
            print *, "[PHYSICS SOLVER] Created:"
            print *, "  Mesh cells: ", mesh%ncells
            print *, "  Domain total cells: ", solver%domain%ntcells
        end if
    end function create_physics_solver
	

	function get_physical_solution(this) result(u_physical)
		class(physics_solver), intent(in) :: this
		real(wp), allocatable :: u_physical(:)
		
		integer :: i, n_physical
		
		n_physical = this%mesh%ncells
		allocate(u_physical(n_physical))
		
		do i = 1, n_physical
			u_physical(i) = this%solution%u(this%domain%ist + i - 1)
		end do
	end function get_physical_solution	
    
    ! ==================== 创建物理组件 ====================
    
    subroutine create_physics_components(this)
        class(physics_solver), intent(inout) :: this
        
        ! 检查是否启用物理模块
        if (.not. this%config%enable_physics) then
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Physics module disabled"
            end if
            return
        end if
        
        ! 创建物理方程
        select case (trim(this%config%equation_type))
        case ("linear_advection")
            allocate(linear_convection_eq :: this%equation)
            select type(eq => this%equation)
            type is(linear_convection_eq)
                eq = create_linear_convection_eq(wave_speed=this%config%wave_speed)
                if (this%config%verbose) then
                    print *, "[PHYSICS SOLVER] Created linear convection equation"
                    print *, "  Wave speed: ", eq%wave_speed
                end if
            end select
            
        case default
            if (this%config%verbose) then
                print *, "[WARNING] Unknown equation type: ", trim(this%config%equation_type)
                print *, "          Using linear convection as default"
            end if
            
            allocate(linear_convection_eq :: this%equation)
            select type(eq => this%equation)
            type is(linear_convection_eq)
                eq = create_linear_convection_eq(wave_speed=this%config%wave_speed)
            end select
        end select
        
        ! 创建物理问题
        select case (trim(this%config%problem_type))
        case ("linear_advection")
            allocate(linear_convection_prob :: this%problem)
            select type(prob => this%problem)
            type is(linear_convection_prob)
                prob = create_linear_convection_prob( &
                    wave_speed=this%config%wave_speed, &
                    domain_length=this%config%domain_length, &
                    ic_type=this%config%ic_type, &
                    boundary_type=this%config%boundary_type)
                if (this%config%verbose) then
                    print *, "[PHYSICS SOLVER] Created linear convection problem"
                    print *, "  IC type: ", trim(prob%ic_type)
                end if
            end select
            
        case default
            if (this%config%verbose) then
                print *, "[WARNING] Unknown problem type: ", trim(this%config%problem_type)
                print *, "          Using linear convection as default"
            end if
            
            allocate(linear_convection_prob :: this%problem)
            select type(prob => this%problem)
            type is(linear_convection_prob)
                prob = create_linear_convection_prob( &
                    wave_speed=this%config%wave_speed, &
                    domain_length=this%config%domain_length, &
                    ic_type=this%config%ic_type, &
                    boundary_type=this%config%boundary_type)
            end select
        end select
        
        this%physics_initialized = .true.
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Physics components created"
        end if
    end subroutine create_physics_components
    
    ! ==================== 创建数值组件 ====================
    
    subroutine create_numerical_components(this)
        class(physics_solver), intent(inout) :: this
        integer :: status
        
        ! 创建重构器
        this%reconstructor = create_reconstructor(this%config, status)
        if (status /= 0) then
            print *, "[ERROR] Failed to create reconstructor"
            this%state = SOLVER_ERROR
            this%error_message = "Failed to create reconstructor"
            return
        end if
        
        ! 创建通量计算器
        this%flux_calculator = create_flux_calculator(this%config, status)
        if (status /= 0) then
            print *, "[ERROR] Failed to create flux calculator"
            this%state = SOLVER_ERROR
            this%error_message = "Failed to create flux calculator"
            return
        end if
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Numerical components created"
        end if
    end subroutine create_numerical_components
    
    ! ==================== 初始化 ====================
    
    subroutine physics_solver_initialize(this)
        class(physics_solver), intent(inout) :: this
        
        if (this%state == SOLVER_INITIALIZED) then
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Already initialized"
            end if
            return
        end if
        
        ! 如果启用了物理模块，应用初始条件
        if (this%physics_initialized .and. allocated(this%problem)) then
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Applying initial condition"
            end if
            
            select type(prob => this%problem)
            type is(linear_convection_prob)
                ! 获取网格单元中心坐标
                call prob%initial_condition(this%mesh%xcc, &
                    this%solution%u(this%domain%ist:this%domain%ied-1))
                
                if (this%config%verbose) then
                    print *, "[PHYSICS SOLVER] Initial condition applied"
                end if
            end select
        else
            ! 简化的初始化：阶跃函数
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Using simplified initialization"
            end if
            
            ! 在 [0.5, 1.0] 区域内设为 2.0，其他区域为 1.0
            where (this%mesh%xcc >= 0.5_wp .and. this%mesh%xcc <= 1.0_wp)
                this%solution%u(this%domain%ist:this%domain%ied-1) = 2.0_wp
            elsewhere
                this%solution%u(this%domain%ist:this%domain%ied-1) = 1.0_wp
            end where
        end if
        
        ! 同步旧场
        call this%solution%update_old_field()
        
        ! 更新状态
        this%state = SOLVER_INITIALIZED
        this%current_time = 0.0_wp
        this%current_step = 0
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Initialized at t = ", this%current_time
        end if
    end subroutine physics_solver_initialize
    
    ! ==================== 时间步进 ====================
    
    subroutine physics_solver_step(this, dt)
        class(physics_solver), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        integer :: i
        real(wp) :: u_val, f_val
        
        if (this%config%verbose .and. mod(this%current_step, 100) == 0) then
            print *, "[PHYSICS SOLVER] Step ", this%current_step + 1, &
                     " dt = ", dt, " t = ", this%current_time
        end if
        
        ! 更新旧场
        call this%solution%update_old_field()
        
        ! 简化的数值方法
        do i = this%domain%ist, this%domain%ied - 1
            u_val = this%solution%un(i)  ! 使用旧值
            
            ! 简单的线性对流：u_t + a*u_x = 0
            ! 使用一阶迎风格式
            this%solution%u(i) = u_val - dt * this%config%wave_speed * &
                (u_val - this%solution%un(i-1)) / this%mesh%dx
        end do
        
        ! 更新时间
        this%current_time = this%current_time + dt
        this%current_step = this%current_step + 1
        
        ! 每100步输出一次进度
        if (this%config%verbose .and. mod(this%current_step, 100) == 0) then
            print *, "[PHYSICS SOLVER] Step ", this%current_step, &
                     " completed, t = ", this%current_time
        end if
    end subroutine physics_solver_step
    
    ! ==================== 运行到指定时间 ====================
    
    subroutine physics_solver_run_to_time(this, final_time)
        class(physics_solver), intent(inout) :: this
        real(wp), intent(in) :: final_time
        
        real(wp) :: dt, t_remaining
        integer :: step_count
        
        if (this%state /= SOLVER_INITIALIZED) then
            this%error_message = "Solver not initialized"
            this%state = SOLVER_ERROR
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER ERROR] Not initialized"
            end if
            return
        end if
        
        this%state = SOLVER_RUNNING
        step_count = 0
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Running from t = ", this%current_time, &
                     " to t = ", final_time
        end if
        
        do while (this%current_time < final_time - 1e-12_wp)
            ! 计算时间步长
            t_remaining = final_time - this%current_time
            dt = min(this%config%dt, t_remaining)
            
            ! 执行时间步
            call this%step(dt)
            
            step_count = step_count + 1
            
            ! 每100步输出一次进度
            if (mod(step_count, 100) == 0 .and. this%config%verbose) then
                print *, "[PHYSICS SOLVER] Progress: t = ", this%current_time, &
                         " / ", final_time
            end if
        end do
        
        ! 恢复原始时间步长
        this%config%dt = this%dt_original
        
        ! 更新状态
        this%state = SOLVER_COMPLETED
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Run completed:"
            print *, "  Final time: ", this%current_time
            print *, "  Total steps: ", this%current_step
            print *, "  Final u range: ", minval(this%solution%u), " to ", maxval(this%solution%u)
        end if
    end subroutine physics_solver_run_to_time
    
    ! ==================== 清理 ====================
    
	subroutine physics_solver_cleanup(this)
		class(physics_solver), intent(inout) :: this
		
		integer :: old_state
		real(wp) :: old_time
		integer(ip) :: old_step
		
		! 保存清理前的状态用于调试
		old_state = this%state
		old_time = this%current_time
		old_step = this%current_step
		
		if (this%config%verbose) then
			print *, "[PHYSICS SOLVER] Cleaning up..."
			print *, "  Before cleanup - State: ", old_state
			print *, "  Before cleanup - Time: ", old_time
			print *, "  Before cleanup - Steps: ", old_step
		end if
		
		! 清理物理组件
		if (allocated(this%equation)) then
			deallocate(this%equation)
			if (this%config%verbose) then
				print *, "  Deallocated equation"
			end if
		end if
		
		if (allocated(this%problem)) then
			deallocate(this%problem)
			if (this%config%verbose) then
				print *, "  Deallocated problem"
			end if
		end if
		
		! 清理数值组件
		if (allocated(this%reconstructor)) then
			deallocate(this%reconstructor)
			if (this%config%verbose) then
				print *, "  Deallocated reconstructor"
			end if
		end if
		
		if (allocated(this%flux_calculator)) then
			deallocate(this%flux_calculator)
			if (this%config%verbose) then
				print *, "  Deallocated flux calculator"
			end if
		end if
		
		! 重置状态 - 但不重置解数组（保持分配）
		this%state = SOLVER_UNINITIALIZED
		this%current_time = 0.0_wp
		this%current_step = 0
		this%error_message = ""
		this%physics_initialized = .false.
		
		if (this%config%verbose) then
			print *, "  After cleanup - State: ", this%state
			print *, "  After cleanup - Time: ", this%current_time
			print *, "  After cleanup - Steps: ", this%current_step
			print *, "[PHYSICS SOLVER] Cleanup completed"
		end if
		
	end subroutine physics_solver_cleanup
    
    ! ==================== 状态查询 ====================
    
    function physics_solver_get_state(this) result(state)
        class(physics_solver), intent(in) :: this
        integer :: state
        state = this%state
    end function physics_solver_get_state
    
    function physics_solver_get_error(this) result(error_msg)
        class(physics_solver), intent(in) :: this
        character(len=100) :: error_msg
        error_msg = trim(this%error_message)
    end function physics_solver_get_error
    
    ! ==================== 信息打印 ====================
    
    subroutine physics_solver_print_info(this)
        class(physics_solver), intent(in) :: this
        
        character(len=20) :: state_str
        
        ! 状态字符串
        select case (this%state)
        case (SOLVER_UNINITIALIZED)
            state_str = "Uninitialized"
        case (SOLVER_INITIALIZED)
            state_str = "Initialized"
        case (SOLVER_RUNNING)
            state_str = "Running"
        case (SOLVER_COMPLETED)
            state_str = "Completed"
        case (SOLVER_ERROR)
            state_str = "Error"
        case default
            write(state_str, '(A, I3)') "Unknown ", this%state
        end select
        
        print *, "=== Physics Solver Information ==="
        print *, "State: ", trim(state_str), " (", this%state, ")"
        print *, "Current time: ", this%current_time
        print *, "Current step: ", this%current_step
        print *, "Error message: '", trim(this%error_message), "'"
        
        ! 配置信息
        print *, "Configuration:"
        print *, "  Scheme: ", trim(this%config%recon_scheme)
        print *, "  Order: ", this%config%spatial_order
        print *, "  dt: ", this%config%dt
        print *, "  Final time: ", this%config%final_time
        
        ! 域信息
        print *, "Domain:"
        print *, "  Ghost layers: ", this%domain%nghosts
        print *, "  Physical cells: ", this%domain%ist, " to ", this%domain%ied - 1
        
        ! 物理信息
        print *, "Physics:"
        print *, "  Initialized: ", this%physics_initialized
        print *, "  Equation type: ", trim(this%config%equation_type)
        print *, "  Problem type: ", trim(this%config%problem_type)
        
        ! 解信息
        if (allocated(this%solution%u)) then
            print *, "Solution:"
            print *, "  u size: ", size(this%solution%u)
            print *, "  u physical size: ", this%domain%ied - this%domain%ist
        end if
        
        print *, "==================================="
    end subroutine physics_solver_print_info

end module physics_solver_module