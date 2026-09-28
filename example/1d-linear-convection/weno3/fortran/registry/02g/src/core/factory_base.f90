! src/core/factory_base.f90
module factory_base_module
    use base_modules, only: wp, ip
    use registry_module, only: create_component, has_component
    
    implicit none
    private
    public :: wp, ip, factory_base, factory_create
    
    ! 工厂基类
    type :: factory_base
        character(len=max_name_length) :: category = ""
    contains
        procedure :: create => factory_base_create
        procedure :: get_available => factory_base_get_available
    end type factory_base
    
    ! 便捷函数类型
    abstract interface
        function factory_function_interface(category, name) result(instance)
            import :: wp
            character(len=*), intent(in) :: category, name
            class(*), allocatable :: instance
        end function factory_function_interface
    end interface
    
contains

    ! 创建工厂实例
    function factory_create(category) result(factory)
        character(len=*), intent(in) :: category
        type(factory_base) :: factory
        factory%category = trim(category)
    end function factory_create
    
    ! 工厂创建方法
    function factory_base_create(this, name) result(instance)
        class(factory_base), intent(in) :: this
        character(len=*), intent(in) :: name
        class(*), allocatable :: instance
        
        instance = create_component(this%category, name)
    end function factory_base_create
    
    ! 获取可用组件列表（简化版）
    subroutine factory_base_get_available(this, names, count)
        class(factory_base), intent(in) :: this
        character(len=*), allocatable, intent(out) :: names(:)
        integer(ip), intent(out) :: count
        
        ! 这里需要实现从注册表获取列表的逻辑
        ! 暂时返回空列表
        count = 0
        allocate(character(len=max_name_length) :: names(0))
    end subroutine factory_base_get_available

end module factory_base_module