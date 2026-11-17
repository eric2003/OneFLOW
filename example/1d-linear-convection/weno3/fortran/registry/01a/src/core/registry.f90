! src/core/registry.f90
module registry_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    
    ! 公开接口
    public :: wp, component_info, component_registry_type, component_registry
    public :: register_component, initialize_registry, cleanup_registry
    
    ! 类型定义
    type :: component_info
        character(len=50) :: category = ""
        character(len=50) :: name = ""
        integer :: order = 0
    contains
        procedure :: print => ci_print
    end type component_info
    
    type :: component_registry_type
        private
        type(component_info), allocatable :: components(:)
        integer :: count = 0
        integer :: capacity = 10
        logical :: verbose = .true.
        logical :: initialized = .false.
    contains
        procedure :: register => cr_register
        procedure :: get => cr_get
        procedure :: list_all => cr_list_all
        procedure :: size => cr_size
        procedure :: clear => cr_clear
    end type component_registry_type
    
    ! 全局实例
    type(component_registry_type) :: component_registry
    
    ! 接口 - 确保两个版本都公开
    interface register_component
        module procedure register_component_simple
        module procedure register_component_with_order
    end interface register_component
    
contains
    
    ! ==================== 公共接口实现 ====================
    
    ! 初始化注册表
    subroutine initialize_registry(initial_capacity, verbose)
        integer, optional, intent(in) :: initial_capacity
        logical, optional, intent(in) :: verbose
        
        if (component_registry%initialized) then
            print *, "[WARN] Registry already initialized"
            return
        end if
        
        if (present(initial_capacity)) then
            component_registry%capacity = max(10, initial_capacity)
        end if
        
        if (present(verbose)) then
            component_registry%verbose = verbose
        end if
        
        component_registry%initialized = .true.
        print *, "[INIT] Registry initialized, capacity:", component_registry%capacity
    end subroutine initialize_registry
    
    ! 清理注册表
    subroutine cleanup_registry
        call component_registry%clear()
        component_registry%initialized = .false.
        print *, "[CLEANUP] Registry cleaned up"
    end subroutine cleanup_registry
    
    ! 简单注册（无阶数）
    subroutine register_component_simple(category, name)
        character(len=*), intent(in) :: category, name
        type(component_info) :: info
        
        info%category = trim(category)
        info%name = trim(name)
        info%order = 0
        
        call component_registry%register(category, name, info)
    end subroutine register_component_simple
    
    ! 带阶数的注册
    subroutine register_component_with_order(category, name, order)
        character(len=*), intent(in) :: category, name
        integer, intent(in) :: order
        type(component_info) :: info
        
        info%category = trim(category)
        info%name = trim(name)
        info%order = order
        
        call component_registry%register(category, name, info)
    end subroutine register_component_with_order
    
    ! ==================== 内部方法 ====================
    
    subroutine ci_print(this)
        class(component_info), intent(in) :: this
        if (this%order > 0) then
            print *, "  [", trim(this%category), ".", trim(this%name), &
                     " (order:", this%order, ")]"
        else
            print *, "  [", trim(this%category), ".", trim(this%name), "]"
        end if
    end subroutine ci_print
    
    ! 内部注册实现
    subroutine cr_register(this, category, name, info)
        class(component_registry_type), intent(inout) :: this
        character(len=*), intent(in) :: category, name
        type(component_info), intent(in) :: info
        
        type(component_info), allocatable :: temp(:)
        integer :: i
        
        if (.not. this%initialized) then
            print *, "[ERROR] Registry not initialized, call initialize_registry first"
            return
        end if
        
        ! 检查是否已存在
        do i = 1, this%count
            if (trim(this%components(i)%category) == trim(category) .and. &
                trim(this%components(i)%name) == trim(name)) then
                if (this%verbose) then
                    print *, "[WARN] Overwriting: ", trim(category), ".", trim(name)
                end if
                this%components(i) = info
                return
            end if
        end do
        
        ! 初始化数组
        if (.not. allocated(this%components)) then
            allocate(this%components(this%capacity))
        end if
        
        ! 扩展数组
        if (this%count >= this%capacity) then
            this%capacity = this%capacity * 2
            allocate(temp(this%capacity))
            temp(1:this%count) = this%components(1:this%count)
            call move_alloc(temp, this%components)
            if (this%verbose) then
                print *, "[INFO] Registry expanded to:", this%capacity
            end if
        end if
        
        ! 添加组件
        this%count = this%count + 1
        this%components(this%count) = info
        
        if (this%verbose) then
            print *, "[OK] Registered: ", trim(category), ".", trim(name)
        end if
    end subroutine cr_register
    
    function cr_get(this, category, name) result(info)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        type(component_info) :: info
        
        integer :: i
        
        info%category = ""
        info%name = ""
        info%order = 0
        
        if (.not. this%initialized) then
            if (this%verbose) then
                print *, "[ERROR] Registry not initialized"
            end if
            return
        end if
        
        do i = 1, this%count
            if (trim(this%components(i)%category) == trim(category) .and. &
                trim(this%components(i)%name) == trim(name)) then
                info = this%components(i)
                return
            end if
        end do
        
        if (this%verbose) then
            print *, "[ERROR] Not found: ", trim(category), ".", trim(name)
        end if
    end function cr_get
    
    subroutine cr_list_all(this)
        class(component_registry_type), intent(in) :: this
        integer :: i
        
        if (.not. this%initialized) then
            print *, "[ERROR] Registry not initialized"
            return
        end if
        
        print *, "=== Registry Contents (", this%count, " components) ==="
        if (this%count == 0) then
            print *, "  (empty)"
            return
        end if
        
        do i = 1, this%count
            call this%components(i)%print()
        end do
    end subroutine cr_list_all
    
    integer function cr_size(this)
        class(component_registry_type), intent(in) :: this
        cr_size = this%count
    end function cr_size
    
    subroutine cr_clear(this)
        class(component_registry_type), intent(inout) :: this
        if (allocated(this%components)) then
            deallocate(this%components)
        end if
        this%count = 0
        this%capacity = 10
        this%initialized = .false.
        print *, "[CLEAR] Registry cleared"
    end subroutine cr_clear
    
    ! 工具函数
    function str(i) result(s)
        integer, intent(in) :: i
        character(len=20) :: s
        write(s, '(I0)') i
        s = trim(adjustl(s))
    end function str
    
end module registry_module