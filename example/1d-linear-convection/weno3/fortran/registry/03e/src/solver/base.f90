! src/solver/base.f90
module solver_base_module
    use base_modules, only: wp => wp, ip => ip  ! 重命名以避免冲突
    
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create
    
    implicit none
    private
    
    ! 明确导出列表
    public :: wp, ip                           ! 类型参数
    public :: solver_base, create_solver_base  ! 类型和构造函数
    public :: SOLVER_UNINITIALIZED, SOLVER_INITIALIZED, SOLVER_RUNNING
    public :: SOLVER_COMPLETED, SOLVER_ERROR   ! 状态常量
    
    ! 求解器状态枚举
    integer, parameter :: SOLVER_UNINITIALIZED = 0
    integer, parameter :: SOLVER_INITIALIZED = 1
    integer, parameter :: SOLVER_RUNNING = 2
    integer, parameter :: SOLVER_COMPLETED = 3
    integer, parameter :: SOLVER_ERROR = 4
    
    ! 求解器基类
    type :: solver_base
        ! 基本组件
        type(cfd_config) :: config
        type(mesh_type) :: mesh
        type(domain_type) :: domain
        type(solution_type) :: solution
        
        ! 状态管理
        integer :: state = SOLVER_UNINITIALIZED
        character(len=100) :: error_message = ""
        real(wp) :: current_time = 0.0_wp
        integer(ip) :: current_step = 0
        
        ! 时间控制
        real(wp) :: dt_original = 0.0_wp
    contains
        procedure :: initialize => solver_base_initialize
        procedure :: step => solver_base_step
        procedure :: run_to_time => solver_base_run_to_time
        procedure :: cleanup => solver_base_cleanup
        procedure :: get_state => solver_base_get_state
        procedure :: get_error => solver_base_get_error
        procedure :: print_info => solver_base_print_info
    end type solver_base
    
    ! 构造函数接口
    interface solver_base
        module procedure create_solver_base
    end interface
    
contains

    ! ==================== 构造函数 ====================
    
    function create_solver_base(config, mesh) result(solver)
        type(cfd_config), intent(in) :: config
        type(mesh_type), intent(in) :: mesh
        type(solver_base) :: solver
        
        solver%config = config
        solver%mesh = mesh
        
        ! 创建域
        solver%domain = domain_create(config, mesh)
        
        ! 创建解
        solver%solution = solution_create(solver%domain)
        
        ! 保存原始时间步长
        solver%dt_original = config%dt
        
        if (config%verbose) then
            print *, "[SOLVER] Base solver created"
            print *, "  Mesh cells: ", mesh%ncells
            print *, "  Domain total cells: ", solver%domain%ntcells
        end if
    end function create_solver_base
    
    ! ==================== 初始化 ====================
    
    subroutine solver_base_initialize(this)
        class(solver_base), intent(inout) :: this
        
        if (this%state == SOLVER_INITIALIZED) then
            if (this%config%verbose) then
                print *, "[SOLVER] Already initialized"
            end if
            return
        end if
        
        ! 初始化解（通过配置）
        ! 这里暂时简化，实际需要调用初始条件工厂
        print *, "[INFO] Base solver initialized (simplified)"
        
        ! 更新状态
        this%state = SOLVER_INITIALIZED
        this%current_time = 0.0_wp
        this%current_step = 0
        
        if (this%config%verbose) then
            print *, "[SOLVER] Initialized at t = ", this%current_time
        end if
    end subroutine solver_base_initialize
    
    ! ==================== 单步计算（虚方法） ====================
    
    subroutine solver_base_step(this, dt)
        class(solver_base), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        ! 基类中这只是虚方法，需要在子类中实现
        print *, "[INFO] Base solver step (virtual method)"
        print *, "  dt = ", dt
        print *, "  t = ", this%current_time
        
        ! 更新时间
        this%current_time = this%current_time + dt
        this%current_step = this%current_step + 1
        
        ! 简单模拟：只是更新状态
        if (this%config%verbose) then
            print *, "[SOLVER] Step completed: t = ", this%current_time, &
                     ", step = ", this%current_step
        end if
    end subroutine solver_base_step
    
    ! ==================== 运行到指定时间 ====================
    
	subroutine solver_base_run_to_time(this, final_time)
		class(solver_base), intent(inout) :: this
		real(wp), intent(in) :: final_time
		
		real(wp) :: dt, t_remaining
		integer :: step_count
		
		if (this%state /= SOLVER_INITIALIZED) then
			this%error_message = "Solver not initialized"
			this%state = SOLVER_ERROR
			if (this%config%verbose) then
				print *, "[SOLVER BASE ERROR] Not initialized: ", trim(this%error_message)
			end if
			return
		end if
		
		this%state = SOLVER_RUNNING
		step_count = 0
		
		if (this%config%verbose) then
			print *, "[SOLVER BASE] Running from t = ", this%current_time, &
					 " to t = ", final_time
			print *, "  Time step: ", this%config%dt
		end if
		
		do while (this%current_time < final_time - 1e-12_wp)
			! 计算时间步长
			t_remaining = final_time - this%current_time
			dt = min(this%config%dt, t_remaining)
			
			! 执行时间步
			call this%step(dt)
			
			step_count = step_count + 1
			
			! 每50步输出一次进度
			if (mod(step_count, 50) == 0 .and. this%config%verbose) then
				print *, "[SOLVER BASE] Progress: t = ", this%current_time, &
						 " / ", final_time, " (step ", step_count, ")"
			end if
		end do
		
		! 恢复原始时间步长
		this%config%dt = this%dt_original
		
		! 更新状态
		this%state = SOLVER_COMPLETED
		
		if (this%config%verbose) then
			print *, "[SOLVER BASE] Run completed:"
			print *, "  Final time: ", this%current_time
			print *, "  Total steps: ", this%current_step
			print *, "  State: ", this%state
		end if
	end subroutine solver_base_run_to_time	
    
    ! ==================== 清理 ====================
    
    subroutine solver_base_cleanup(this)
        class(solver_base), intent(inout) :: this
        
        ! 重置状态
        this%state = SOLVER_UNINITIALIZED
        this%current_time = 0.0_wp
        this%current_step = 0
        this%error_message = ""
        
        if (this%config%verbose) then
            print *, "[SOLVER] Cleaned up"
        end if
    end subroutine solver_base_cleanup
    
    ! ==================== 状态查询 ====================
    
    function solver_base_get_state(this) result(state)
        class(solver_base), intent(in) :: this
        integer :: state
        state = this%state
    end function solver_base_get_state
    
    function solver_base_get_error(this) result(error_msg)
        class(solver_base), intent(in) :: this
        character(len=100) :: error_msg
        error_msg = trim(this%error_message)
    end function solver_base_get_error
    
    ! ==================== 信息打印 ====================
    
    subroutine solver_base_print_info(this)
        class(solver_base), intent(in) :: this
        
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
            state_str = "Unknown"
        end select
        
        print *, "=== Solver Information ==="
        print *, "State: ", trim(state_str)
        print *, "Current time: ", this%current_time
        print *, "Current step: ", this%current_step
        print *, "Error message: '", trim(this%error_message), "'"
        
        ! 配置信息
        print *, "Configuration:"
        print *, "  Scheme: ", trim(this%config%recon_scheme)
        print *, "  Order: ", this%config%spatial_order
        print *, "  dt: ", this%config%dt
        
        ! 域信息
        print *, "Domain:"
        print *, "  Ghost layers: ", this%domain%nghosts
        print *, "  Physical cells: ", this%domain%ist, " to ", this%domain%ied - 1
        
        print *, "========================="
    end subroutine solver_base_print_info

end module solver_base_module