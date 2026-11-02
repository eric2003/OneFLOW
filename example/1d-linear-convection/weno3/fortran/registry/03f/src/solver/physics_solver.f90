! src/solver/physics_solver.f90 (修正版)
module physics_solver_module
    use base_modules, only: wp => wp, ip => ip
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create
    
    implicit none
    private
    
    ! 明确导出列表
    public :: wp, ip, physics_solver
    public :: SOLVER_UNINITIALIZED, SOLVER_INITIALIZED, SOLVER_COMPLETED, SOLVER_ERROR
    
    ! 求解器状态枚举
    integer, parameter :: SOLVER_UNINITIALIZED = 0
    integer, parameter :: SOLVER_INITIALIZED = 1
    integer, parameter :: SOLVER_RUNNING = 2
    integer, parameter :: SOLVER_COMPLETED = 3
    integer, parameter :: SOLVER_ERROR = 4
    
    ! 物理求解器类型
    type :: physics_solver
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
        procedure :: initialize => physics_solver_initialize
        procedure :: run_to_time => physics_solver_run_to_time
        procedure :: cleanup => physics_solver_cleanup
        procedure :: get_state => physics_solver_get_state
        procedure, private :: apply_simple_initial_condition  ! 添加这行
    end type physics_solver
    
contains

    ! ==================== 初始化 ====================
    
    subroutine physics_solver_initialize(this)
        class(physics_solver), intent(inout) :: this
        
        if (this%state == SOLVER_INITIALIZED) then
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Already initialized"
            end if
            return
        end if
        
        ! 创建域
        this%domain = domain_create(this%config, this%mesh)
        
        ! 创建解
        this%solution = solution_create(this%domain)
        
        ! 应用简化的初始条件
        call this%apply_simple_initial_condition()
        
        ! 保存原始时间步长
        this%dt_original = this%config%dt
        
        ! 更新状态
        this%state = SOLVER_INITIALIZED
        this%current_time = 0.0_wp
        this%current_step = 0
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Initialized at t = ", this%current_time
        end if
    end subroutine physics_solver_initialize
    
    subroutine apply_simple_initial_condition(this)
        class(physics_solver), intent(inout) :: this
        
        integer :: i, idx
        real(wp) :: x
        
        ! 简化的阶跃函数初始条件
        do i = this%domain%ist, this%domain%ied - 1
            idx = i - this%domain%ist + 1
            x = this%mesh%xcc(idx)
            
            if (x >= 0.5_wp .and. x <= 1.0_wp) then
                this%solution%u(i) = 2.0_wp
            else
                this%solution%u(i) = 1.0_wp
            end if
        end do
        
        ! 同步旧场
        this%solution%un = this%solution%u
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Applied simple initial condition"
        end if
    end subroutine
    
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
            
            ! 简单的时间步进（占位符）
            this%current_time = this%current_time + dt
            this%current_step = this%current_step + 1
            step_count = step_count + 1
            
            ! 每100步输出一次进度
            if (mod(step_count, 100) == 0 .and. this%config%verbose) then
                print *, "[PHYSICS SOLVER] Progress: t = ", this%current_time, &
                         " / ", final_time, " (step ", step_count, ")"
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
        end if
    end subroutine physics_solver_run_to_time
    
    ! ==================== 清理 ====================
    
    subroutine physics_solver_cleanup(this)
        class(physics_solver), intent(inout) :: this
        
        this%state = SOLVER_UNINITIALIZED
        this%current_time = 0.0_wp
        this%current_step = 0
        this%error_message = ""
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Cleaned up"
        end if
    end subroutine physics_solver_cleanup
    
    ! ==================== 状态查询 ====================
    
    function physics_solver_get_state(this) result(state)
        class(physics_solver), intent(in) :: this
        integer :: state
        state = this%state
    end function physics_solver_get_state
    
end module physics_solver_module