! src/core/registry_advanced.f90
module registry_advanced_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use factory_interfaces, only: factory_procedure
    implicit none
    
    private
    public :: wp, component_info, component_registry, register_factory, create_instance
    
    ! 组件信息类型（增强版）
    type :: component_info
        character(len=50) :: category = ""
        character(len=50) :: name = ""
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
        integer :: capacity = 10
        logical :: verbose = .true.
    contains
        procedure :: register => cr_register
        procedure :: get => cr_get
        procedure :: list_all => cr_list_all
        procedure :: create => cr_create
    end type component_registry_type
    
    type(component_registry_type) :: component_registry
    
    ! 接口
    interface register_factory
        module procedure register_factory_simple
        module procedure register_factory_with_order
    end interface register_factory
    
contains
    
    ! ==================== 组件信息方法 ====================
    
    subroutine ci_print(this)
        class(component_info), intent(in) :: this
        character(len=100) :: factory_str
        
        if (this%has_factory) then
            factory_str = " [has factory]"
        else
            factory_str = ""
        end if
        
        if (this%order > 0) then
            print *, "  [", trim(this%category), ".", trim(this%name), &
                     " (order:", this%order, ")", trim(factory_str), "]"
        else
            print *, "  [", trim(this%category), ".", trim(this%name), &
                     trim(factory_str), "]"
        end if
    end subroutine ci_print
    
    subroutine ci_create(this, instance)
        class(component_info), intent(in) :: this
        class(*), allocatable, intent(out) :: instance
        
        if (.not. this%has_factory) then
            error stop "[ERROR] Component has no factory"
        end if
        
        call this%factory(instance)
    end subroutine ci_create
    
    ! ==================== 注册表方法 ====================
    
    ! 注册工厂（简单版）
    subroutine register_factory_simple(category, name, factory_proc)
        character(len=*), intent(in) :: category, name
        procedure(factory_procedure) :: factory_proc
        
        type(component_info) :: info
        
        info%category = trim(category)
        info%name = trim(name)
        info%order = 0
        info%factory => factory_proc
        info%has_factory = .true.
        
        call component_registry%register(info)
    end subroutine register_factory_simple
    
    ! 注册工厂（带阶数）
    subroutine register_factory_with_order(category, name, factory_proc, order)
        character(len=*), intent(in) :: category, name
        procedure(factory_procedure) :: factory_proc
        integer, intent(in) :: order
        
        type(component_info) :: info
        
        info%category = trim(category)
        info%name = trim(name)
        info%order = order
        info%factory => factory_proc
        info%has_factory = .true.
        
        call component_registry%register(info)
    end subroutine register_factory_with_order
    
    ! 内部注册实现
    subroutine cr_register(this, info)
        class(component_registry_type), intent(inout) :: this
        type(component_info), intent(in) :: info
        
        type(component_info), allocatable :: temp(:)
        integer :: i
        
        ! 检查是否已存在
        do i = 1, this%count
            if (trim(this%components(i)%category) == trim(info%category) .and. &
                trim(this%components(i)%name) == trim(info%name)) then
                if (this%verbose) then
                    print *, "[WARN] Overwriting: ", &
                             trim(info%category), ".", trim(info%name)
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
        end if
        
        ! 添加组件
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
        
        ! 初始化为空
        info%category = ""
        info%name = ""
        info%order = 0
        info%factory => null()
        info%has_factory = .false.
        
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
    
    ! 创建实例
    subroutine cr_create(this, category, name, instance)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        class(*), allocatable, intent(out) :: instance
        
        type(component_info) :: info
        
        info = this%get(category, name)
        
        if (len_trim(info%category) > 0 .and. info%has_factory) then
            call info%create(instance)
        else
            error stop "[ERROR] Cannot create instance: factory not available"
        end if
    end subroutine cr_create
    
    subroutine cr_list_all(this)
        class(component_registry_type), intent(in) :: this
        integer :: i
        
        print *, "=== Registry Contents (", this%count, " components) ==="
        if (this%count == 0) then
            print *, "  (empty)"
            return
        end if
        
        do i = 1, this%count
            call this%components(i)%print()
        end do
    end subroutine cr_list_all
    
end module registry_advanced_module