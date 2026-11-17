! src/numerics/reconstructor/factory.f90
module reconstructor_factory_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use registry_module, only: create_component, has_component
    use reconstructor_base_module, only: reconstructor_base
    implicit none
    
    private
    public :: wp, reconstructor_factory_create
    
contains
    
    subroutine reconstructor_factory_create(scheme_name, reconstructor)
        character(len=*), intent(in) :: scheme_name
        class(reconstructor_base), allocatable, intent(out) :: reconstructor
        
        class(*), allocatable :: instance
        character(len=20) :: name_lower
        integer :: i
        
        ! 转换为小写
        name_lower = scheme_name
        do i = 1, len_trim(name_lower)
            if (name_lower(i:i) >= 'A' .and. name_lower(i:i) <= 'Z') then
                name_lower(i:i) = char(ichar(name_lower(i:i)) + 32)
            end if
        end do
        
        ! 检查是否已注册
        if (.not. has_component("reconstructor", trim(name_lower))) then
            print *, "[ERROR] 重构器未注册: ", trim(scheme_name)
            return
        end if
        
        ! 创建实例
        call create_component("reconstructor", trim(name_lower), instance)
        
        ! 类型转换
        select type (inst => instance)
        type is (reconstructor_base)
            allocate(reconstructor, source=inst)
        class default
            error stop "[ERROR] 创建的重构器类型错误"
        end select
    end subroutine reconstructor_factory_create
    
end module reconstructor_factory_module