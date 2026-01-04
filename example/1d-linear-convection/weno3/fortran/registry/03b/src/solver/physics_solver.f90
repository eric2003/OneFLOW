! src/solver/physics_solver.f90 (修复版)
module physics_solver_module
    use base_modules, only: wp => wp, ip => ip
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use solver_base_module, only: solver_base, create_solver_base
    use solver_base_module, only: SOLVER_UNINITIALIZED, SOLVER_INITIALIZED, &
                                 SOLVER_RUNNING, SOLVER_COMPLETED, SOLVER_ERROR
    
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
    
    ! 基于物理的求解器（扩展自solver_base）
    type, extends(solver_base) :: physics_solver
        ! 物理组件
        class(physics_equation), allocatable :: equation
        class(physics_problem), allocatable :: problem
        
        ! 数值组件
        class(reconstructor_base), allocatable :: reconstructor
        class(flux_calculator_base), allocatable :: flux_calculator
        
        ! 状态信息
        logical :: physics_initialized = .false.
    contains
        procedure :: initialize => physics_solver_initialize
        procedure :: step => physics_solver_step
        procedure :: cleanup => physics_solver_cleanup
        procedure :: print_info => physics_solver_print_info
        procedure :: create_physics_components
        procedure :: create_numerical_components
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
        
        ! 先调用父类构造函数
        solver%solver_base = create_solver_base(config, mesh)
        
        ! 创建物理组件
        call solver%create_physics_components()
        
        ! 创建数值组件
        call solver%create_numerical_components()
        
        if (config%verbose) then
            print *, "[PHYSICS SOLVER] Created with physics support"
        end if
    end function create_physics_solver
    
    ! ==================== 创建物理组件 ====================
    
    subroutine create_physics_components(this)
        class(physics_solver), intent(inout) :: this
        
        ! 检查是否启用物理模块
        if (.not. this%config%enable_physics) then
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Physics module disabled, skipping physics components"
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
            print *, "[WARNING] Unknown equation type: ", trim(this%config%equation_type)
            print *, "          Using linear convection as default"
            
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
                    print *, "  Domain length: ", prob%domain_length
                end if
            end select
            
        case default
            print *, "[WARNING] Unknown problem type: ", trim(this%config%problem_type)
            print *, "          Using linear convection as default"
            
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
            print *, "[PHYSICS SOLVER] Physics components created successfully"
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
            print *, "[PHYSICS SOLVER] Numerical components created successfully"
        end if
    end subroutine create_numerical_components
    
    ! ==================== 初始化 ====================
    
    subroutine physics_solver_initialize(this)
        class(physics_solver), intent(inout) :: this
        
        ! 调用父类初始化
        call this%solver_base%initialize()
        
        ! 如果启用了物理模块，应用初始条件
        if (this%physics_initialized .and. allocated(this%problem)) then
            if (this%config%verbose) then
                print *, "[PHYSICS SOLVER] Applying initial condition from physics problem"
            end if
            
            select type(prob => this%problem)
            type is(linear_convection_prob)
                ! 获取网格单元中心坐标
                call prob%initial_condition(this%mesh%xcc, this%solution%u(this%domain%ist:this%domain%ied-1))
                
                if (this%config%verbose) then
                    print *, "[PHYSICS SOLVER] Initial condition applied"
                    print *, "  u range: ", minval(this%solution%u), " to ", maxval(this%solution%u)
                end if
            end select
        end if
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Initialized with physics support"
        end if
    end subroutine physics_solver_initialize
    
    ! ==================== 时间步进（简化实现） ====================
    
    subroutine physics_solver_step(this, dt)
        class(physics_solver), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        integer :: i, j
        real(wp) :: u_val, f_val, residual
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Step with dt = ", dt
        end if
        
        ! 简化的数值方法：直接使用方程计算通量
        if (allocated(this%equation) .and. allocated(this%reconstructor) .and. &
            allocated(this%flux_calculator)) then
            
            ! 这里应该是完整的数值方法：
            ! 1. 边界条件
            ! 2. 重构
            ! 3. 通量计算
            ! 4. 残差计算
            ! 5. 时间积分
            
            ! 简化的占位实现：只是模拟计算过程
            do i = this%domain%ist, this%domain%ied - 1
                j = i - this%domain%ist + 1
                u_val = this%solution%u(i)
                
                ! 使用方程计算通量（简化）
                select type(eq => this%equation)
                type is(linear_convection_eq)
                    f_val = eq%flux(u_val)
                end select
                
                ! 简化的更新：只是演示
                this%solution%u(i) = u_val - dt * f_val / this%mesh%dx
            end do
            
            ! 更新时间
            this%current_time = this%current_time + dt
            this%current_step = this%current_step + 1
            
            if (this%config%verbose .and. mod(this%current_step, 10) == 0) then
                print *, "[PHYSICS SOLVER] Step ", this%current_step, &
                         " completed, t = ", this%current_time
            end if
            
        else
            ! 如果没有完整的组件，回退到基类方法
            call this%solver_base%step(dt)
        end if
    end subroutine physics_solver_step
    
    ! ==================== 清理 ====================
    
    subroutine physics_solver_cleanup(this)
        class(physics_solver), intent(inout) :: this
        
        ! 清理物理组件
        if (allocated(this%equation)) deallocate(this%equation)
        if (allocated(this%problem)) deallocate(this%problem)
        
        ! 清理数值组件
        if (allocated(this%reconstructor)) deallocate(this%reconstructor)
        if (allocated(this%flux_calculator)) deallocate(this%flux_calculator)
        
        ! 调用父类清理
        call this%solver_base%cleanup()
        
        this%physics_initialized = .false.
        
        if (this%config%verbose) then
            print *, "[PHYSICS SOLVER] Cleaned up physics components"
        end if
    end subroutine physics_solver_cleanup
    
    ! ==================== 信息打印 ====================
    
    subroutine physics_solver_print_info(this)
        class(physics_solver), intent(in) :: this
        
        ! 调用父类信息打印
        call this%solver_base%print_info()
        
        ! 添加物理信息
        print *, "=== Physics Information ==="
        print *, "Physics initialized: ", this%physics_initialized
        
        if (allocated(this%equation)) then
            print *, "Equation: allocated"
            select type(eq => this%equation)
            type is(linear_convection_eq)
                print *, "  Type: Linear Convection"
                print *, "  Wave speed: ", eq%wave_speed
            class default
                print *, "  Type: Unknown"
            end select
        else
            print *, "Equation: not allocated"
        end if
        
        if (allocated(this%problem)) then
            print *, "Problem: allocated"
            select type(prob => this%problem)
            type is(linear_convection_prob)
                print *, "  Type: Linear Convection Problem"
                print *, "  IC type: ", trim(prob%ic_type)
                print *, "  Domain length: ", prob%domain_length
            class default
                print *, "  Type: Unknown"
            end select
        else
            print *, "Problem: not allocated"
        end if
        
        print *, "Numerical components:"
        print *, "  Reconstructor: ", allocated(this%reconstructor)
        print *, "  Flux calculator: ", allocated(this%flux_calculator)
        print *, "=========================="
    end subroutine physics_solver_print_info

end module physics_solver_module