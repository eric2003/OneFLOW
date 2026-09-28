! src/core/registry.f90
module registry_module
    use, intrinsic :: iso_fortran_env, only: real64
    use factory_interfaces, only: factory_procedure
    implicit none

    private

    ! ==================== PUBLIC INTERFACE ====================
    ! 原有接口保持不变
    public :: real64, component_info, component_registry
    public :: register_component_simple, initialize_registry, cleanup_registry
    public :: has_component, get_available_components
    public :: registry_is_initialized, registry_get_size
    public :: list_all  ! ← 恢复这个接口
    
    ! 新增工厂相关接口
    public :: factory_ptr_type, component_registry_entry
    public :: register_component_with_factory
    public :: create_component_from_registry
    public :: has_factory_component
    public :: list_factory_components

    ! ==================== TYPE DEFINITIONS ====================
    
    ! 原有类型定义
    type :: component_info
        character(len=32) :: category = ""
        character(len=32) :: name = ""
        integer :: order = 0
    contains
        procedure :: print => ci_print
    end type component_info
    
    ! 新增：工厂指针类型
    type :: factory_ptr_type
        procedure(factory_procedure), pointer, nopass :: ptr => null()
    end type factory_ptr_type
    
    ! 新增：带工厂的注册条目
    type :: component_registry_entry
        character(len=32) :: category = ""
        character(len=32) :: name = ""
        integer :: order = 0
        type(factory_ptr_type) :: factory
        logical :: has_factory = .false.
    contains
        procedure :: print_with_factory => entry_print
        procedure :: can_create => entry_has_factory
    end type component_registry_entry

    ! 原有注册表类型（扩展以支持工厂）
    type :: component_registry_type
        private
        ! 简单注册组件（向后兼容）
        type(component_info), allocatable :: simple_components(:)
        integer :: simple_count = 0
        
        ! 带工厂的组件（新增）
        type(component_registry_entry), allocatable :: factory_components(:)
        integer :: factory_count = 0
        
        integer :: capacity = 100
        logical :: verbose = .true.
        logical :: initialized = .false.
    contains
        ! 简单注册方法
        procedure :: register => cr_register           ! ← 保持原名
        procedure :: get => cr_get                     ! ← 保持原名
        procedure :: list_all => cr_list_all           ! ← 保持原名（重要！）
        procedure :: list_simple => cr_list_simple     ! ← 新增别名
        
        ! 工厂注册方法
        procedure :: register_with_factory => cr_register_with_factory
        procedure :: get_with_factory => cr_get_with_factory
        procedure :: list_with_factory => cr_list_with_factory
        procedure :: create_from_factory => cr_create_from_factory
        
        ! 通用方法
        procedure :: clear => cr_clear
        procedure :: size => cr_size
        procedure :: total_size => cr_total_size
        procedure :: is_initialized => cr_is_initialized
    end type component_registry_type

    ! Global registry instance
    type(component_registry_type), save :: component_registry

contains

    ! ==================== PUBLIC API ====================

    ! Initialize registry (保持不变)
    subroutine initialize_registry(initial_capacity, verbose)
        integer, optional, intent(in) :: initial_capacity
        logical, optional, intent(in) :: verbose

        if (component_registry%initialized) then
            if (component_registry%verbose) then
                print *, "[INFO] Registry already initialized"
            end if
            return
        end if

        if (present(initial_capacity)) then
            component_registry%capacity = max(10, initial_capacity)
        end if

        if (present(verbose)) then
            component_registry%verbose = verbose
        end if

        ! Allocate arrays for both types
        allocate(component_registry%simple_components(component_registry%capacity))
        allocate(component_registry%factory_components(component_registry%capacity))

        component_registry%initialized = .true.
        component_registry%simple_count = 0
        component_registry%factory_count = 0

        if (component_registry%verbose) then
            print *, "[INIT] Registry initialized, capacity:", component_registry%capacity
            print *, "       Supports both simple and factory-based registration"
        end if
    end subroutine initialize_registry

    ! Cleanup registry (保持不变)
    subroutine cleanup_registry
        call component_registry%clear()
        if (component_registry%verbose) then
            print *, "[CLEANUP] Registry cleaned up"
        end if
    end subroutine cleanup_registry

    ! ==================== SIMPLE REGISTRATION (向后兼容) ====================

    ! Simple registration (no factory) - 保持不变
    subroutine register_component_simple(category, name)
        character(len=*), intent(in) :: category, name

        type(component_info) :: info

        info%category = to_lower(trim(adjustl(category)))
        info%name = to_lower(trim(adjustl(name)))
        info%order = 0

        call component_registry%register(info)
    end subroutine register_component_simple

    ! Check if component exists (简单注册)
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

    ! Get available components in a category (简单注册)
    subroutine get_available_components(category, names, orders)
        character(len=*), intent(in) :: category
        character(len=:), allocatable, intent(out), optional :: names(:)
        integer, allocatable, intent(out), optional :: orders(:)

        character(len=32) :: cat_lower
        integer :: i, count, idx
        type(component_info) :: info

        cat_lower = to_lower(trim(adjustl(category)))

        ! Count components in this category
        count = 0
        do i = 1, component_registry%simple_count
            if (component_registry%simple_components(i)%category == cat_lower) then
                count = count + 1
            end if
        end do

        ! Allocate arrays if requested
        if (present(names)) then
            allocate(character(len=32) :: names(count))
        end if

        if (present(orders)) then
            allocate(orders(count))
        end if

        ! Fill arrays
        idx = 1
        do i = 1, component_registry%simple_count
            if (component_registry%simple_components(i)%category == cat_lower) then
                info = component_registry%simple_components(i)
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

    ! ==================== FACTORY REGISTRATION (新增) ====================

    ! Register component with factory function
    subroutine register_component_with_factory(category, name, factory_func, order)
        character(len=*), intent(in) :: category, name
        procedure(factory_procedure) :: factory_func
        integer, optional, intent(in) :: order
        
        type(component_registry_entry) :: entry
        character(len=32) :: cat_lower, name_lower
        
        cat_lower = to_lower(trim(adjustl(category)))
        name_lower = to_lower(trim(adjustl(name)))
        
        entry%category = cat_lower
        entry%name = name_lower
        
        if (present(order)) then
            entry%order = order
        else
            entry%order = 0
        end if
        
        entry%factory%ptr => factory_func
        entry%has_factory = .true.
        
        call component_registry%register_with_factory(entry)
    end subroutine register_component_with_factory

    ! Create component from registry using factory
    function create_component_from_registry(category, name) result(instance)
        character(len=*), intent(in) :: category, name
        class(*), allocatable :: instance
        
        character(len=32) :: cat_lower, name_lower
        integer :: status
        
        cat_lower = to_lower(trim(adjustl(category)))
        name_lower = to_lower(trim(adjustl(name)))
        
        call component_registry%create_from_factory(cat_lower, name_lower, instance, status)
        
        if (status /= 0) then
            if (component_registry%verbose) then
                print *, "[ERROR] Failed to create component: ", trim(category), ".", trim(name)
            end if
        end if
    end function create_component_from_registry

    ! Check if component has factory
    function has_factory_component(category, name) result(found)
        character(len=*), intent(in) :: category, name
        logical :: found
        
        character(len=32) :: cat_lower, name_lower
        
        cat_lower = to_lower(trim(adjustl(category)))
        name_lower = to_lower(trim(adjustl(name)))
        
        found = .false.
        
        if (.not. component_registry%initialized) return
        
        ! 简化的查找逻辑（实际应调用内部方法）
        found = component_registry%factory_count > 0  ! 占位实现
    end function has_factory_component

    ! List all factory components
    subroutine list_factory_components(category)
        character(len=*), optional, intent(in) :: category
        
        if (present(category)) then
            call component_registry%list_with_factory(category)
        else
            call component_registry%list_with_factory("")
        end if
    end subroutine list_factory_components

    ! Public function to list all simple components (向后兼容)
    subroutine list_all()
        call component_registry%list_all()
    end subroutine list_all

    ! ==================== PUBLIC UTILITY FUNCTIONS ====================

    ! Public function to check if registry is initialized
    function registry_is_initialized() result(is_initialized)
        logical :: is_initialized
        is_initialized = component_registry%is_initialized()
    end function registry_is_initialized

    ! Public function to get total registry size
    function registry_get_size() result(size_val)
        integer :: size_val
        size_val = component_registry%total_size()
    end function registry_get_size

    ! ==================== COMPONENT INFO METHODS ====================

    subroutine ci_print(this)
        class(component_info), intent(in) :: this

        if (this%order > 0) then
            print *, "  [", trim(this%category), ".", trim(this%name), &
                     " (order:", this%order, ")]"
        else
            print *, "  [", trim(this%category), ".", trim(this%name), "]"
        end if
    end subroutine ci_print

    subroutine entry_print(this)
        class(component_registry_entry), intent(in) :: this
        
        if (this%has_factory) then
            if (this%order > 0) then
                print *, "  [FACTORY] ", trim(this%category), ".", trim(this%name), &
                         " (order:", this%order, ")"
            else
                print *, "  [FACTORY] ", trim(this%category), ".", trim(this%name)
            end if
        else
            if (this%order > 0) then
                print *, "  [SIMPLE] ", trim(this%category), ".", trim(this%name), &
                         " (order:", this%order, ")"
            else
                print *, "  [SIMPLE] ", trim(this%category), ".", trim(this%name)
            end if
        end if
    end subroutine entry_print

    logical function entry_has_factory(this)
        class(component_registry_entry), intent(in) :: this
        entry_has_factory = this%has_factory
    end function entry_has_factory

    ! ==================== REGISTRY INTERNAL METHODS ====================

    ! ---------- Simple Registration Methods ----------
    
    subroutine cr_register(this, info)
        class(component_registry_type), intent(inout) :: this
        type(component_info), intent(in) :: info

        type(component_info), allocatable :: temp(:)
        integer :: i

        if (.not. this%initialized) then
            error stop "[ERROR] Registry not initialized, call initialize_registry first"
        end if

        ! Check if already exists
        do i = 1, this%simple_count
            if (this%simple_components(i)%category == info%category .and. &
                this%simple_components(i)%name == info%name) then
                if (this%verbose) then
                    print *, "[WARN] Overwriting simple component: ", &
                             trim(info%category), ".", trim(info%name)
                end if
                this%simple_components(i) = info
                return
            end if
        end do

        ! Expand array if needed
        if (this%simple_count >= this%capacity) then
            this%capacity = this%capacity * 2
            allocate(temp(this%capacity))
            temp(1:this%simple_count) = this%simple_components(1:this%simple_count)
            call move_alloc(temp, this%simple_components)

            if (this%verbose) then
                print *, "[INFO] Simple registry expanded to capacity:", this%capacity
            end if
        end if

        ! Add component
        this%simple_count = this%simple_count + 1
        this%simple_components(this%simple_count) = info

        if (this%verbose) then
            print *, "[OK] Registered simple: ", trim(info%category), ".", trim(info%name)
        end if
    end subroutine cr_register

    function cr_get(this, category, name) result(info)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        type(component_info) :: info

        integer :: i

        ! Initialize return value as empty
        info%category = ""
        info%name = ""
        info%order = 0

        if (.not. this%initialized) then
            return
        end if

        do i = 1, this%simple_count
            if (this%simple_components(i)%category == category .and. &
                this%simple_components(i)%name == name) then
                info = this%simple_components(i)
                return
            end if
        end do
    end function cr_get

    subroutine cr_list_all(this)
        class(component_registry_type), intent(in) :: this
        integer :: i

        if (.not. this%initialized) then
            print *, "[INFO] Registry not initialized"
            return
        end if

        if (this%simple_count == 0) then
            print *, "[INFO] No simple components registered"
            return
        end if

        print *, "=== Simple Registry Contents (", this%simple_count, " components) ==="
        do i = 1, this%simple_count
            call this%simple_components(i)%print()
        end do
        print *, "=========================================="
    end subroutine cr_list_all

    subroutine cr_list_simple(this)
        class(component_registry_type), intent(in) :: this
        call this%list_all()
    end subroutine cr_list_simple

    ! ---------- Factory Registration Methods ----------
    
    subroutine cr_register_with_factory(this, entry)
        class(component_registry_type), intent(inout) :: this
        type(component_registry_entry), intent(in) :: entry

        type(component_registry_entry), allocatable :: temp(:)
        integer :: i

        if (.not. this%initialized) then
            error stop "[ERROR] Registry not initialized, call initialize_registry first"
        end if

        ! Check if already exists
        do i = 1, this%factory_count
            if (this%factory_components(i)%category == entry%category .and. &
                this%factory_components(i)%name == entry%name) then
                if (this%verbose) then
                    print *, "[WARN] Overwriting factory component: ", &
                             trim(entry%category), ".", trim(entry%name)
                end if
                this%factory_components(i) = entry
                return
            end if
        end do

        ! Expand array if needed
        if (this%factory_count >= this%capacity) then
            this%capacity = this%capacity * 2
            allocate(temp(this%capacity))
            temp(1:this%factory_count) = this%factory_components(1:this%factory_count)
            call move_alloc(temp, this%factory_components)

            if (this%verbose) then
                print *, "[INFO] Factory registry expanded to capacity:", this%capacity
            end if
        end if

        ! Add component
        this%factory_count = this%factory_count + 1
        this%factory_components(this%factory_count) = entry

        if (this%verbose) then
            print *, "[FACTORY] Registered with factory: ", &
                     trim(entry%category), ".", trim(entry%name)
        end if
    end subroutine cr_register_with_factory

    function cr_get_with_factory(this, category, name) result(entry)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        type(component_registry_entry) :: entry

        integer :: i

        ! Initialize return value as empty
        entry%category = ""
        entry%name = ""
        entry%order = 0
        entry%factory%ptr => null()
        entry%has_factory = .false.

        if (.not. this%initialized) then
            return
        end if

        do i = 1, this%factory_count
            if (this%factory_components(i)%category == category .and. &
                this%factory_components(i)%name == name) then
                entry = this%factory_components(i)
                return
            end if
        end do
    end function cr_get_with_factory

    subroutine cr_list_with_factory(this, category)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category
        
        character(len=32) :: cat_lower
        integer :: i, count
        
        if (.not. this%initialized) then
            print *, "[INFO] Registry not initialized"
            return
        end if
        
        cat_lower = to_lower(trim(adjustl(category)))
        
        if (this%factory_count == 0) then
            print *, "[INFO] No factory components registered"
            return
        end if
        
        ! Count components in this category
        count = 0
        do i = 1, this%factory_count
            if (len_trim(cat_lower) == 0 .or. &
                this%factory_components(i)%category == cat_lower) then
                count = count + 1
            end if
        end do
        
        if (count == 0) then
            if (len_trim(cat_lower) > 0) then
                print *, "[INFO] No factory components in category: ", trim(cat_lower)
            else
                print *, "[INFO] No factory components registered"
            end if
            return
        end if
        
        print *, "=== Factory Registry Contents (", count, " components) ==="
        do i = 1, this%factory_count
            if (len_trim(cat_lower) == 0 .or. &
                this%factory_components(i)%category == cat_lower) then
                call this%factory_components(i)%print_with_factory()
            end if
        end do
        print *, "=========================================="
    end subroutine cr_list_with_factory

    subroutine cr_create_from_factory(this, category, name, instance, status)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        class(*), allocatable, intent(out) :: instance
        integer, intent(out) :: status
        
        type(component_registry_entry) :: entry
        integer :: i
        
        status = -1  ! Default: not found
        
        if (.not. this%initialized) then
            if (this%verbose) then
                print *, "[ERROR] Registry not initialized"
            end if
            return
        end if
        
        ! Find the entry
        do i = 1, this%factory_count
            if (this%factory_components(i)%category == category .and. &
                this%factory_components(i)%name == name) then
                entry = this%factory_components(i)
                
                if (entry%has_factory .and. associated(entry%factory%ptr)) then
                    ! Call the factory function
                    call entry%factory%ptr(instance)
                    status = 0  ! Success
                    
                    if (this%verbose) then
                        print *, "[FACTORY] Created component: ", &
                                 trim(category), ".", trim(name)
                    end if
                else
                    status = -2  ! No factory function
                    if (this%verbose) then
                        print *, "[ERROR] No factory function for: ", &
                                 trim(category), ".", trim(name)
                    end if
                end if
                return
            end if
        end do
        
        ! If we get here, component not found
        if (this%verbose) then
            print *, "[ERROR] Factory component not found: ", &
                     trim(category), ".", trim(name)
        end if
    end subroutine cr_create_from_factory

    ! ---------- Common Methods ----------
    
    subroutine cr_clear(this)
        class(component_registry_type), intent(inout) :: this

        if (allocated(this%simple_components)) then
            deallocate(this%simple_components)
        end if
        
        if (allocated(this%factory_components)) then
            deallocate(this%factory_components)
        end if

        this%simple_count = 0
        this%factory_count = 0
        this%capacity = 100
        this%initialized = .false.
    end subroutine cr_clear

    integer function cr_total_size(this)
        class(component_registry_type), intent(in) :: this
        cr_total_size = this%simple_count + this%factory_count
    end function cr_total_size
    
    integer function cr_size(this)
        class(component_registry_type), intent(in) :: this
        cr_size = this%simple_count  ! Backward compatibility
    end function cr_size

    logical function cr_is_initialized(this)
        class(component_registry_type), intent(in) :: this
        cr_is_initialized = this%initialized
    end function cr_is_initialized

    ! ==================== UTILITY FUNCTIONS ====================

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