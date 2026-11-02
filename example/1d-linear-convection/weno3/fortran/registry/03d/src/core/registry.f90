! src/core/registry.f90 (更新版)
module registry_module
    use base_modules, only: wp, ip, max_name_len, component_info
    
    implicit none
    private
    
    ! 明确公开所有需要的接口
    public :: wp, ip                            ! 类型参数
    public :: component_info                    ! 类型
    public :: registry_init, registry_cleanup   ! 初始化/清理
    public :: register_component_simple         ! 注册组件
    public :: has_component_simple              ! 检查组件
    public :: list_components                   ! 列出组件
    public :: registry_is_initialized           ! 检查初始化状态
    public :: registry_get_size                 ! 获取大小
    public :: initialize_default_components     ! 新增：初始化默认组件
    
    ! 全局注册表
    type :: component_registry
        type(component_info), allocatable :: components(:)
        integer(ip) :: count = 0
        integer(ip) :: capacity = 100
        logical :: initialized = .false.
        logical :: verbose = .true.
        logical :: default_components_added = .false.  ! 新增：标记是否已添加默认组件
    end type component_registry
    
    type(component_registry) :: registry
    
contains

    ! ==================== 公共API ====================
    
    subroutine registry_init(verbose)
        logical, optional, intent(in) :: verbose
        
        if (registry%initialized) then
            if (registry%verbose) then
                print *, "[REGISTRY] Already initialized"
            end if
            return
        end if
        
        if (present(verbose)) then
            registry%verbose = verbose
        end if
        
        allocate(registry%components(registry%capacity))
        registry%initialized = .true.
        
        if (registry%verbose) then
            print *, "[REGISTRY] Initialized with capacity:", registry%capacity
        end if
    end subroutine registry_init
    
    subroutine registry_cleanup()
        if (allocated(registry%components)) then
            deallocate(registry%components)
        end if
        registry%initialized = .false.
        registry%count = 0
        registry%default_components_added = .false.  ! 重置标记
        
        if (registry%verbose) then
            print *, "[REGISTRY] Cleaned up"
        end if
    end subroutine registry_cleanup
    
    ! 新增：初始化默认组件
    subroutine initialize_default_components()
        if (.not. registry%initialized) then
            call registry_init()
        end if
        
        if (registry%default_components_added) then
            if (registry%verbose) then
                print *, "[REGISTRY] Default components already added"
            end if
            return
        end if
        
        ! 注册重构器
        call register_component_simple("reconstructor", "eno", order=3)
        call register_component_simple("reconstructor", "weno3", order=3)
        call register_component_simple("reconstructor", "weno5", order=5)
        
        ! 注册通量计算器
        call register_component_simple("flux", "rusanov")
        call register_component_simple("flux", "engquist-osher")
        
        ! 注册边界条件
        call register_component_simple("boundary", "periodic")
        call register_component_simple("boundary", "dirichlet")
        call register_component_simple("boundary", "neumann")
        
        ! 注册时间积分器
        call register_component_simple("integrator", "rk1", order=1)
        call register_component_simple("integrator", "rk2", order=2)
        call register_component_simple("integrator", "rk3", order=3)
        
        ! 注册方程
        call register_component_simple("equation", "linear_advection")
        
        ! 注册问题
        call register_component_simple("problem", "linear_advection")
        
        ! 注册初始条件
        call register_component_simple("initial_condition", "step")
        call register_component_simple("initial_condition", "sin")
        call register_component_simple("initial_condition", "gaussian")
        
        registry%default_components_added = .true.
        
        if (registry%verbose) then
            print *, "[REGISTRY] Default components registered"
            print *, "[REGISTRY] Total components:", registry%count
        end if
    end subroutine initialize_default_components
    
    subroutine register_component_simple(category, name, order)
        character(len=*), intent(in) :: category, name
        integer(ip), optional, intent(in) :: order
        
        integer(ip) :: i
        type(component_info) :: info
        
        if (.not. registry%initialized) then
            call registry_init()
        end if
        
        ! 检查是否已存在
        do i = 1, registry%count
            if (trim(registry%components(i)%category) == trim(category) .and. &
                trim(registry%components(i)%name) == trim(name)) then
                if (registry%verbose) then
                    print *, "[WARN] Overwriting component: ", trim(category), ".", trim(name)
                end if
                
                ! 更新
                if (present(order)) then
                    registry%components(i)%order = order
                else
                    registry%components(i)%order = 0
                end if
                return
            end if
        end do
        
        ! 扩展数组
        if (registry%count >= registry%capacity) then
            call expand_registry()
        end if
        
        ! 添加新组件
        registry%count = registry%count + 1
        
        info%category = trim(category)
        info%name = trim(name)
        info%order = 0
        if (present(order)) then
            info%order = order
        end if
        
        registry%components(registry%count) = info
        
        if (registry%verbose) then
            print *, "[OK] Registered simple: ", trim(category), ".", trim(name)
        end if
    end subroutine register_component_simple
    
    logical function has_component_simple(category, name)
        character(len=*), intent(in) :: category, name
        
        integer(ip) :: i
        
        has_component_simple = .false.
        
        if (.not. registry%initialized) return
        
        do i = 1, registry%count
            if (trim(registry%components(i)%category) == trim(category) .and. &
                trim(registry%components(i)%name) == trim(name)) then
                has_component_simple = .true.
                return
            end if
        end do
    end function has_component_simple
    
    subroutine list_components(category)
        character(len=*), optional, intent(in) :: category
        
        integer(ip) :: i, count
        
        if (.not. registry%initialized) then
            print *, "[INFO] Registry not initialized"
            return
        end if
        
        if (registry%count == 0) then
            print *, "[INFO] No components registered"
            return
        end if
        
        count = 0
        print *, "=== Registry Contents ==="
        do i = 1, registry%count
            if (.not. present(category) .or. &
                trim(registry%components(i)%category) == trim(category)) then
                call print_component_info(registry%components(i))
                count = count + 1
            end if
        end do
        
        print *, "Total:", count, "components"
        print *, "=========================="
    end subroutine list_components
    
    ! ==================== 新增函数 ====================
    
    logical function registry_is_initialized()
        ! 检查注册表是否已初始化
        registry_is_initialized = registry%initialized
    end function registry_is_initialized
    
    integer(ip) function registry_get_size()
        ! 获取注册表中的组件数量
        registry_get_size = registry%count
    end function registry_get_size
    
    ! ==================== 内部辅助函数 ====================
    
    subroutine expand_registry()
        type(component_info), allocatable :: temp(:)
        
        registry%capacity = registry%capacity * 2
        allocate(temp(registry%capacity))
        temp(1:registry%count) = registry%components(1:registry%count)
        call move_alloc(temp, registry%components)
        
        if (registry%verbose) then
            print *, "[INFO] Registry expanded to capacity:", registry%capacity
        end if
    end subroutine expand_registry
    
    subroutine print_component_info(info)
        type(component_info), intent(in) :: info
        
        if (info%order > 0) then
            print *, "  [", trim(info%category), ".", trim(info%name), &
                     " (order:", info%order, ")]"
        else
            print *, "  [", trim(info%category), ".", trim(info%name), "]"
        end if
    end subroutine print_component_info

end module registry_module