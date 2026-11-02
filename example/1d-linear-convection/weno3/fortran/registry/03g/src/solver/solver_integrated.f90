! src/solver/solver_integrated.f90 (完全修正版)
module solver_integrated_module
    use base_modules, only: wp, ip
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create
    use boundary_base_module, only: boundary_condition
    use periodic_boundary_module, only: periodic_boundary
    
    implicit none
    private
    
    ! 求解器状态常量
    integer, parameter, public :: SOLVER_READY = 0
    integer, parameter, public :: SOLVER_INITIALIZED = 1
    integer, parameter, public :: SOLVER_RUNNING = 2
    integer, parameter, public :: SOLVER_COMPLETED = 3
    integer, parameter, public :: SOLVER_ERROR = -1
    
    type, public :: integrated_solver
        ! 基本组件
        type(cfd_config) :: config
        type(mesh_type) :: mesh
        type(domain_type) :: domain
        type(solution_type) :: solution
        
        ! 组件实例
        class(boundary_condition), allocatable :: bc
        
        ! 状态
        integer :: state = SOLVER_READY
        real(wp) :: current_time = 0.0_wp
        integer(ip) :: current_step = 0
        character(len=100) :: error_msg = ""
        
        ! 数据模式
        logical :: use_real_data = .true.
        
    contains
        procedure :: initialize => solver_initialize
        procedure :: run_to_time => solver_run_to_time
        procedure :: cleanup => solver_cleanup
        procedure :: get_state => solver_get_state
        procedure :: enable_real_data => solver_enable_real_data
        procedure, private :: apply_initial_condition
        procedure, private :: apply_boundary_conditions
        procedure, private :: simple_time_step
        procedure, private :: calculate_dt
        procedure, private :: setup_boundary_condition  ! 修改方法名
    end type integrated_solver
    
contains

    ! ==================== 公共接口 ====================
    
    subroutine solver_initialize(this)
        class(integrated_solver), intent(inout) :: this
        
        if (this%state /= SOLVER_READY) then
            this%error_msg = "Solver already initialized"
            this%state = SOLVER_ERROR
            if (this%config%verbose) then
                print *, "[ERROR] ", trim(this%error_msg)
            end if
            return
        end if
        
        ! 创建域
        this%domain = domain_create(this%config, this%mesh)
        
        ! 创建解
        this%solution = solution_create(this%domain)
        
        ! 设置边界条件
        call this%setup_boundary_condition()  ! 修改调用
        
        ! 应用初始条件
        call this%apply_initial_condition()
        
        ! 初始化状态
        this%state = SOLVER_INITIALIZED
        this%current_time = 0.0_wp
        this%current_step = 0
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Initialized successfully"
            print *, "  Domain cells (with ghosts): ", this%domain%ntcells
            print *, "  Ghost layers: ", this%domain%nghosts
            if (this%use_real_data) then
                print *, "  Data mode: Real"
            else
                print *, "  Data mode: Simple"
            end if
        end if
    end subroutine solver_initialize
    
    ! ==================== 边界条件设置 ====================
    
    subroutine setup_boundary_condition(this)
        class(integrated_solver), intent(inout) :: this
        
        ! 根据配置创建边界条件
        select case (trim(this%config%boundary_type))
        case ("periodic")
            ! 创建周期性边界条件
            if (.not. allocated(this%bc)) then
                allocate(periodic_boundary :: this%bc)
            end if
            
            select type(bc => this%bc)
            type is (periodic_boundary)
                bc%name = "periodic"
            end select
            
        case default
            ! 默认使用周期性边界
            if (this%config%verbose) then
                print *, "[WARNING] Using periodic boundary as default"
            end if
            
            if (.not. allocated(this%bc)) then
                allocate(periodic_boundary :: this%bc)
            end if
            
            select type(bc => this%bc)
            type is (periodic_boundary)
                bc%name = "periodic"
            end select
        end select
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Boundary condition created: ", &
                     trim(this%bc%get_name())
        end if
    end subroutine setup_boundary_condition
    
    subroutine solver_run_to_time(this, final_time)
        class(integrated_solver), intent(inout) :: this
        real(wp), intent(in) :: final_time
        
        real(wp) :: dt, t_remaining
        integer :: step_count
        real(wp) :: original_dt
        
        if (this%state /= SOLVER_INITIALIZED) then
            this%error_msg = "Solver not initialized"
            this%state = SOLVER_ERROR
            if (this%config%verbose) then
                print *, "[ERROR] ", trim(this%error_msg)
            end if
            return
        end if
        
        ! 保存原始时间步长
        original_dt = this%config%dt
        
        this%state = SOLVER_RUNNING
        step_count = 0
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Starting time integration"
            print *, "  Initial time: ", this%current_time
            print *, "  Final time: ", final_time
            print *, "  Initial dt: ", this%config%dt
        end if
        
        do while (this%current_time < final_time - 1e-12_wp)
            ! 计算时间步长
            call this%calculate_dt()
            
            ! 确保不超过最终时间
            t_remaining = final_time - this%current_time
            dt = min(this%config%dt, t_remaining)
            
            ! 应用边界条件
            call this%apply_boundary_conditions()
            
            ! 执行时间步
            call this%simple_time_step(dt)
            
            ! 更新时间
            this%current_time = this%current_time + dt
            this%current_step = this%current_step + 1
            step_count = step_count + 1
            
            ! 进度输出
            if (mod(step_count, 50) == 0 .and. this%config%verbose) then
                print *, "  Progress: t = ", this%current_time, &
                         " / ", final_time, " (step ", step_count, ")"
            end if
        end do
        
        ! 恢复原始时间步长
        this%config%dt = original_dt
        
        ! 更新状态
        this%state = SOLVER_COMPLETED
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Time integration completed"
            print *, "  Final time: ", this%current_time
            print *, "  Total steps: ", this%current_step
            if (allocated(this%solution%u)) then
                print *, "  Solution range: [", minval(this%solution%u), ", ", &
                        maxval(this%solution%u), "]"
            end if
        end if
    end subroutine solver_run_to_time
    
    subroutine solver_cleanup(this)
        class(integrated_solver), intent(inout) :: this
        
        ! 清理分配的组件
        if (allocated(this%bc)) then
            deallocate(this%bc)
        end if
        
        ! 重置状态
        this%state = SOLVER_READY
        this%current_time = 0.0_wp
        this%current_step = 0
        this%error_msg = ""
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Cleaned up"
        end if
    end subroutine solver_cleanup
    
    function solver_get_state(this) result(state)
        class(integrated_solver), intent(in) :: this
        integer :: state
        state = this%state
    end function solver_get_state
    
    subroutine solver_enable_real_data(this, use_real)
        class(integrated_solver), intent(inout) :: this
        logical, intent(in) :: use_real
        this%use_real_data = use_real
        
        if (this%config%verbose) then
            if (use_real) then
                print *, "[INTEGRATED SOLVER] Data mode set to: Real"
            else
                print *, "[INTEGRATED SOLVER] Data mode set to: Simple"
            end if
        end if
    end subroutine solver_enable_real_data
    
    ! ==================== 私有方法 ====================
    
    subroutine apply_initial_condition(this)
        class(integrated_solver), intent(inout) :: this
        
        integer :: i, idx
        real(wp) :: x
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Applying initial condition: ", &
                     trim(this%config%ic_type)
        end if
        
        ! 简化的初始条件应用
        do i = this%domain%ist, this%domain%ied - 1
            idx = i - this%domain%ist + 1
            x = this%mesh%xcc(idx)
            
            select case (trim(this%config%ic_type))
            case ("step")
                ! 阶跃函数
                if (x >= 0.5_wp .and. x <= 1.0_wp) then
                    this%solution%u(i) = 2.0_wp
                else
                    this%solution%u(i) = 1.0_wp
                end if
                
            case ("sin", "sine")
                ! 正弦波
                this%solution%u(i) = sin(2.0_wp * 3.141592653589793_wp * x / &
                                       this%config%domain_length)
                
            case ("gaussian")
                ! 高斯脉冲
                this%solution%u(i) = exp(-((x - 0.5_wp) / 0.1_wp)**2)
                
            case default
                ! 默认阶跃函数
                if (x >= 0.5_wp .and. x <= 1.0_wp) then
                    this%solution%u(i) = 2.0_wp
                else
                    this%solution%u(i) = 1.0_wp
                end if
            end select
        end do
        
        ! 同步旧场
        this%solution%un = this%solution%u
        
        if (this%config%verbose) then
            print *, "[INTEGRATED SOLVER] Initial condition applied"
            if (allocated(this%solution%u)) then
                print *, "  Min value: ", minval(this%solution%u)
                print *, "  Max value: ", maxval(this%solution%u)
            end if
        end if
    end subroutine apply_initial_condition
    
    subroutine apply_boundary_conditions(this)
        class(integrated_solver), intent(inout) :: this
        
        if (.not. allocated(this%bc)) then
            if (this%config%verbose) then
                print *, "[WARNING] No boundary condition allocated"
            end if
            return
        end if
        
        select type(bc => this%bc)
        type is (periodic_boundary)
            ! 应用周期性边界条件
            call bc%apply(this%solution%u, &
                          this%domain%nghosts, &
                          this%domain%ist, &
                          this%domain%ied - 1)
            
        class default
            ! 对于其他边界条件类型
            if (this%config%verbose) then
                print *, "[WARNING] Boundary condition type not fully implemented"
            end if
        end select
    end subroutine apply_boundary_conditions
    
    subroutine simple_time_step(this, dt)
        class(integrated_solver), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        integer :: i
        real(wp) :: dx, cfl
        
        ! 简化的时间步进（一阶迎风格式）
        dx = this%mesh%dx
        cfl = this%config%wave_speed * dt / dx
        
        if (cfl > 1.0_wp .and. this%config%verbose) then
            print *, "[WARNING] CFL = ", cfl, " > 1.0"
        end if
        
        ! 保存旧解
        this%solution%un = this%solution%u
        
        ! 一阶迎风格式（只更新内部点）
        if (this%domain%ist + 1 <= this%domain%ied - 2) then
            do i = this%domain%ist + 1, this%domain%ied - 2
                this%solution%u(i) = this%solution%un(i) - &
                    cfl * (this%solution%un(i) - this%solution%un(i-1))
            end do
        end if
        
        ! 调试输出
        if (mod(this%current_step, 100) == 0 .and. this%config%verbose) then
            print *, "  [TIME STEP] Step: ", this%current_step, &
                     ", t = ", this%current_time, &
                     ", CFL = ", cfl, &
                     ", dt = ", dt
        end if
    end subroutine simple_time_step
    
    subroutine calculate_dt(this)
        class(integrated_solver), intent(inout) :: this
        
        real(wp) :: cfl, dx
        
        dx = this%mesh%dx
        
        if (this%use_real_data) then
            ! 真实计算使用CFL条件
            cfl = 0.8_wp  ! CFL数
            this%config%dt = cfl * dx / abs(this%config%wave_speed)
        else
            ! 简单数据使用固定时间步长
            this%config%dt = 0.0025_wp
        end if
        
        if (this%config%verbose .and. this%current_step == 0) then
            print *, "[INTEGRATED SOLVER] Calculated dt = ", this%config%dt
        end if
    end subroutine calculate_dt

end module solver_integrated_module