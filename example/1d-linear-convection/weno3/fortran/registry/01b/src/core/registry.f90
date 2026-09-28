! src/core/registry.f90
module registry_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    
    ! ==================== 公开接口 ====================
    public :: wp, component_info, component_registry
    public :: register_component, create_component
    public :: initialize_registry, cleanup_registry
    public :: has_component, get_available_components
    
    ! ==================== 类型定义 ====================
    
    ! 工厂过程接口
    abstract interface
        subroutine factory_procedure(instance)
            import :: wp
            class(*), allocatable, intent(out) :: instance
        end subroutine factory_procedure
    end interface
    
    ! 组件信息类型
    type :: component_info
        character(len=32) :: category = ""
        character(len=32) :: name = ""
        integer :: order = 0
        procedure(factory_procedure), pointer, nopass :: factory => null()
        logical :: has_factory = .false.
    contains
        procedure :: print => ci_print
        procedure :: create => ci_create
    end type component_info
    
    ! 注册表类型
    type :: component_registry_type
        private
        type(component_info), allocatable :: components(:)
        integer :: count = 0
        integer :: capacity = 100
        logical :: verbose = .true.
        logical :: initialized = .false.
    contains
        procedure, private :: register_info => cr_register_info
        procedure :: get => cr_get
        procedure :: list_all => cr_list_all
        procedure :: clear => cr_clear
        procedure :: size => cr_size
        procedure, private :: expand => cr_expand
    end type component_registry_type
    
    ! 全局注册表实例
    type(component_registry_type), save :: component_registry
    
    ! 接口重载
    interface register_component
        module procedure register_component_simple
        module procedure register_component_with_factory
        module procedure register_component_full
    end interface register_component
    
contains
    
    ! ==================== 公共API ====================
    
    ! 初始化注册表
    subroutine initialize_registry(initial_capacity, verbose)
        integer, optional, intent(in) :: initial_capacity
        logical, optional, intent(in) :: verbose
        
        if (component_registry%initialized) then
            if (component_registry%verbose) then
                print *, "[INFO] 注册表已初始化"
            end if
            return
        end if
        
        if (present(initial_capacity)) then
            component_registry%capacity = max(10, initial_capacity)
        end if
        
        if (present(verbose)) then
            component_registry%verbose = verbose
        end if
        
        ! 分配数组
        allocate(component_registry%components(component_registry%capacity))
        
        component_registry%initialized = .true.
        component_registry%count = 0
        
        if (component_registry%verbose) then
            print *, "[INIT] 注册表初始化完成，容量:", component_registry%capacity
        end if
    end subroutine initialize_registry
    
    ! 清理注册表
    subroutine cleanup_registry
        call component_registry%clear()
        if (component_registry%verbose) then
            print *, "[CLEANUP] 注册表清理完成"
        end if
    end subroutine cleanup_registry
    
    ! 简单注册（只有名称）
    subroutine register_component_simple(category, name)
        character(len=*), intent(in) :: category, name
        type(component_info) :: info
        
        info%category = to_lower(trim(adjustl(category)))
        info%name = to_lower(trim(adjustl(name)))
        info%order = 0
        info%factory => null()
        info%has_factory = .false.
        
        call component_registry%register_info(info)
    end subroutine register_component_simple
    
    ! 带阶数的注册
    subroutine register_component_with_factory(category, name, factory_proc, order)
        character(len=*), intent(in) :: category, name
        procedure(factory_procedure) :: factory_proc
        integer, optional, intent(in) :: order
        
        type(component_info) :: info
        
        info%category = to_lower(trim(adjustl(category)))
        info%name = to_lower(trim(adjustl(name)))
        
        if (present(order)) then
            info%order = order
        else
            info%order = 0
        end if
        
        info%factory => factory_proc
        info%has_factory = .true.
        
        call component_registry%register_info(info)
    end subroutine register_component_with_factory
    
    ! 完整注册
    subroutine register_component_full(category, name, order, has_factory)
        character(len=*), intent(in) :: category, name
        integer, intent(in) :: order
        logical, intent(in) :: has_factory
        
        type(component_info) :: info
        
        info%category = to_lower(trim(adjustl(category)))
        info%name = to_lower(trim(adjustl(name)))
        info%order = order
        info%factory => null()
        info%has_factory = has_factory
        
        call component_registry%register_info(info)
    end subroutine register_component_full
    
    ! 创建组件实例
    subroutine create_component(category, name, instance)
        character(len=*), intent(in) :: category, name
        class(*), allocatable, intent(out) :: instance
        
        type(component_info) :: info
        character(len=32) :: cat_lower, name_lower
        
        cat_lower = to_lower(trim(adjustl(category)))
        name_lower = to_lower(trim(adjustl(name)))
        
        info = component_registry%get(cat_lower, name_lower)
        
        if (len_trim(info%category) == 0) then
            error stop "[ERROR] 组件未找到: " // trim(cat_lower) // "." // trim(name_lower)
        end if
        
        if (.not. info%has_factory) then
            error stop "[ERROR] 组件没有工厂函数: " // trim(cat_lower) // "." // trim(name_lower)
        end if
        
        if (.not. associated(info%factory)) then
            error stop "[ERROR] 工厂函数未关联: " // trim(cat_lower) // "." // trim(name_lower)
        end if
        
        call info%create(instance)
    end subroutine create_component
    
    ! 检查组件是否存在
    function has_component(category, name) result(found)
        character(len=*), intent(in) :: category, name
        logical :: found
        
        type(component_info) :: info
        character(len=32) :: cat_lower, name_lower
        
        cat_lower = to_lower(trim(adjustl(category)))
        name_lower = to_lower(trim(adjustl(name)))
        
        info = component_registry%get(cat_lower, name_lower)
        found = (len_trim(info%category) > 0)
    end function has_component
    
    ! 获取某类别下的所有组件
    subroutine get_available_components(category, names, orders)
        character(len=*), intent(in) :: category
        character(len=:), allocatable, intent(out), optional :: names(:)
        integer, allocatable, intent(out), optional :: orders(:)
        
        character(len=32) :: cat_lower
        integer :: i, count, idx
        type(component_info) :: info
        
        cat_lower = to_lower(trim(adjustl(category)))
        
        ! 先计算数量
        count = 0
        do i = 1, component_registry%count
            if (component_registry%components(i)%category == cat_lower) then
                count = count + 1
            end if
        end do
        
        ! 分配数组
        if (present(names)) then
            allocate(character(len=32) :: names(count))
        end if
        
        if (present(orders)) then
            allocate(orders(count))
        end if
        
        ! 填充数组
        idx = 1
        do i = 1, component_registry%count
            if (component_registry%components(i)%category == cat_lower) then
                info = component_registry%components(i)
                if (present(names)) then
                    names(idx) = info%name
                end if
                if (present(orders)) then
                    orders(idx) = info%order
                end if
                idx = idx + 1
            end if
        end do
    end subroutine get_available_components
    
    ! ==================== 组件信息方法 ====================
    
    subroutine ci_print(this)
        class(component_info), intent(in) :: this
        
        if (this%order > 0) then
            if (this%has_factory) then
                print *, "  [", trim(this%category), ".", trim(this%name), &
                         " (order:", this%order, ", has factory)]"
            else
                print *, "  [", trim(this%category), ".", trim(this%name), &
                         " (order:", this%order, ")]"
            end if
        else
            if (this%has_factory) then
                print *, "  [", trim(this%category), ".", trim(this%name), &
                         " (has factory)]"
            else
                print *, "  [", trim(this%category), ".", trim(this%name), "]"
            end if
        end if
    end subroutine ci_print
    
    subroutine ci_create(this, instance)
        class(component_info), intent(in) :: this
        class(*), allocatable, intent(out) :: instance
        
        if (.not. this%has_factory) then
            error stop "[ERROR] 组件没有工厂函数"
        end if
        
        if (.not. associated(this%factory)) then
            error stop "[ERROR] 工厂函数未关联"
        end if
        
        call this%factory(instance)
    end subroutine ci_create
    
    ! ==================== 注册表内部方法 ====================
    
    subroutine cr_register_info(this, info)
        class(component_registry_type), intent(inout) :: this
        type(component_info), intent(in) :: info
        
        integer :: i
        
        if (.not. this%initialized) then
            error stop "[ERROR] 注册表未初始化，请先调用 initialize_registry"
        end if
        
        ! 检查是否已存在
        do i = 1, this%count
            if (this%components(i)%category == info%category .and. &
                this%components(i)%name == info%name) then
                if (this%verbose) then
                    print *, "[WARN] 覆盖注册: ", &
                             trim(info%category), ".", trim(info%name)
                end if
                this%components(i) = info
                return
            end if
        end do
        
        ! 检查是否需要扩展
        if (this%count >= this%capacity) then
            call this%expand()
        end if
        
        ! 添加组件
        this%count = this%count + 1
        this%components(this%count) = info
        
        if (this%verbose) then
            print *, "[OK] 已注册: ", trim(info%category), ".", trim(info%name)
        end if
    end subroutine cr_register_info
    
    function cr_get(this, category, name) result(info)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        type(component_info) :: info
        
        integer :: i
        
        ! 初始化返回值为空
        info%category = ""
        info%name = ""
        info%order = 0
        info%factory => null()
        info%has_factory = .false.
        
        if (.not. this%initialized) then
            return
        end if
        
        do i = 1, this%count
            if (this%components(i)%category == category .and. &
                this%components(i)%name == name) then
                info = this%components(i)
                return
            end if
        end do
    end function cr_get
    
    subroutine cr_list_all(this)
        class(component_registry_type), intent(in) :: this
        integer :: i
        
        if (.not. this%initialized) then
            print *, "[INFO] 注册表未初始化"
            return
        end if
        
        print *, "=== 注册表内容 (", this%count, " 个组件) ==="
        
        if (this%count == 0) then
            print *, "  （空）"
            return
        end if
        
        ! 按类别分组显示
        do i = 1, this%count
            call this%components(i)%print()
        end do
        
        print *, "======================================"
    end subroutine cr_list_all
    
    subroutine cr_clear(this)
        class(component_registry_type), intent(inout) :: this
        
        if (allocated(this%components)) then
            deallocate(this%components)
        end if
        
        this%count = 0
        this%capacity = 100
        this%initialized = .false.
    end subroutine cr_clear
    
    integer function cr_size(this)
        class(component_registry_type), intent(in) :: this
        cr_size = this%count
    end function cr_size
    
    subroutine cr_expand(this)
        class(component_registry_type), intent(inout) :: this
        type(component_info), allocatable :: temp(:)
        
        this%capacity = this%capacity * 2
        allocate(temp(this%capacity))
        temp(1:this%count) = this%components(1:this%count)
        call move_alloc(temp, this%components)
        
        if (this%verbose) then
            print *, "[INFO] 注册表扩展至容量:", this%capacity
        end if
    end subroutine cr_expand
    
    ! ==================== 工具函数 ====================
    
    function to_lower(str) result(lower_str)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: lower_str
        integer :: i
        
        do i = 1, len(str)
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
                lower_str(i:i) = char(ichar(str(i:i)) + 32)
            else
                lower_str(i:i) = str(i:i)
            end if
        end do
    end function to_lower
    
end module registry_module