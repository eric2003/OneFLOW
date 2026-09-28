! src/core/registry.f90
module registry_module
    use, intrinsic :: iso_fortran_env, only: real64
    use factory_interfaces, only: factory_procedure
    implicit none
    
    private
    
    ! Public interface
    public :: real64, component_info, component_registry
    public :: register_component_simple, initialize_registry, cleanup_registry
    public :: has_component, get_available_components
    public :: registry_is_initialized, registry_get_size  ! 添加公共访问方法
    
    ! Type definitions (simplified, no factory for now)
    type :: component_info
        character(len=32) :: category = ""
        character(len=32) :: name = ""
        integer :: order = 0
    contains
        procedure :: print => ci_print
    end type component_info
    
    type :: component_registry_type
        private
        type(component_info), allocatable :: components(:)
        integer :: count = 0
        integer :: capacity = 100
        logical :: verbose = .true.
        logical :: initialized = .false.
    contains
        procedure :: register => cr_register
        procedure :: get => cr_get
        procedure :: list_all => cr_list_all
        procedure :: clear => cr_clear
        procedure :: size => cr_size
        procedure :: is_initialized => cr_is_initialized  ! 内部方法
    end type component_registry_type
    
    ! Global registry instance
    type(component_registry_type), save :: component_registry
    
contains
    
    ! ==================== PUBLIC API ====================
    
    ! Initialize registry
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
        
        ! Allocate array
        allocate(component_registry%components(component_registry%capacity))
        
        component_registry%initialized = .true.
        component_registry%count = 0
        
        if (component_registry%verbose) then
            print *, "[INIT] Registry initialized, capacity:", component_registry%capacity
        end if
    end subroutine initialize_registry
    
    ! Cleanup registry
    subroutine cleanup_registry
        call component_registry%clear()
        if (component_registry%verbose) then
            print *, "[CLEANUP] Registry cleaned up"
        end if
    end subroutine cleanup_registry
    
    ! Simple registration (no factory)
    subroutine register_component_simple(category, name)
        character(len=*), intent(in) :: category, name
        
        type(component_info) :: info
        
        info%category = to_lower(trim(adjustl(category)))
        info%name = to_lower(trim(adjustl(name)))
        info%order = 0
        
        call component_registry%register(info)
    end subroutine register_component_simple
    
    ! Check if component exists
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
    
    ! Get available components in a category
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
        do i = 1, component_registry%count
            if (component_registry%components(i)%category == cat_lower) then
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
    
    ! Public function to check if registry is initialized
    function registry_is_initialized() result(is_initialized)
        logical :: is_initialized
        is_initialized = component_registry%is_initialized()
    end function registry_is_initialized
    
    ! Public function to get registry size
    function registry_get_size() result(size_val)
        integer :: size_val
        size_val = component_registry%size()
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
    
    ! ==================== REGISTRY INTERNAL METHODS ====================
    
    subroutine cr_register(this, info)
        class(component_registry_type), intent(inout) :: this
        type(component_info), intent(in) :: info
        
        type(component_info), allocatable :: temp(:)
        integer :: i
        
        if (.not. this%initialized) then
            error stop "[ERROR] Registry not initialized, call initialize_registry first"
        end if
        
        ! Check if already exists
        do i = 1, this%count
            if (this%components(i)%category == info%category .and. &
                this%components(i)%name == info%name) then
                if (this%verbose) then
                    print *, "[WARN] Overwriting: ", &
                             trim(info%category), ".", trim(info%name)
                end if
                this%components(i) = info
                return
            end if
        end do
        
        ! Expand array if needed
        if (this%count >= this%capacity) then
            this%capacity = this%capacity * 2
            allocate(temp(this%capacity))
            temp(1:this%count) = this%components(1:this%count)
            call move_alloc(temp, this%components)
            
            if (this%verbose) then
                print *, "[INFO] Registry expanded to capacity:", this%capacity
            end if
        end if
        
        ! Add component
        this%count = this%count + 1
        this%components(this%count) = info
        
        if (this%verbose) then
            print *, "[OK] Registered: ", trim(info%category), ".", trim(info%name)
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
            print *, "[INFO] Registry not initialized"
            return
        end if
        
        print *, "=== Registry Contents (", this%count, " components) ==="
        
        if (this%count == 0) then
            print *, "  (empty)"
            return
        end if
        
        ! Show components grouped by category
        do i = 1, this%count
            call this%components(i)%print()
        end do
        
        print *, "==========================================="
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