! src/core/factory_interfaces.f90
module factory_interfaces
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    public :: wp
    
    ! 抽象工厂接口
    type, abstract :: base_factory
    contains
        procedure(create_interface), deferred :: create
    end type base_factory
    
    abstract interface
        subroutine create_interface(this, instance)
            import :: base_factory
            class(base_factory), intent(in) :: this
            class(*), allocatable, intent(out) :: instance
        end subroutine create_interface
    end interface
    
end module factory_interfaces