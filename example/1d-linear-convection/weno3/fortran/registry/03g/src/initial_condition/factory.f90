! src/initial_condition/factory.f90
module ic_factory_module
    use base_modules, only: wp, ip
    use ic_base_module, only: initial_condition
    use step_ic_module, only: step_function_ic
    implicit none
    private
    
    public :: create_initial_condition, initialize_ic_factory
    
contains

    subroutine initialize_ic_factory()
        print *, "[IC FACTORY] Initial conditions placeholder"
    end subroutine
    
    subroutine create_initial_condition(ic_type, ic_instance)
        character(len=*), intent(in) :: ic_type
        class(initial_condition), allocatable, intent(out) :: ic_instance
        
        ! 暂时只实现步函数
        allocate(step_function_ic :: ic_instance)
        
        select type(ic => ic_instance)
        type is (step_function_ic)
            ! 使用默认构造函数
            ic%name = "step"
        end select
    end subroutine
    
end module ic_factory_module