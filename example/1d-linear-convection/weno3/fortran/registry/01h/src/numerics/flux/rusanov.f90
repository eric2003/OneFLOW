! src/numerics/flux/rusanov.f90
module rusanov_flux_module
    use, intrinsic :: iso_fortran_env, only: real64
    use flux_base_module, only: flux_calculator_base
    implicit none
    
    private
    public :: real64, rusanov_flux
    
    type, extends(flux_calculator_base) :: rusanov_flux
        real(real64) :: wave_speed_default = 1.0_real64
    contains
        procedure :: info => rusanov_info
    end type rusanov_flux
    
    ! 添加构造函数接口
    interface rusanov_flux
        module procedure create_rusanov_flux
    end interface
    
contains
    
    ! 构造函数
    type(rusanov_flux) function create_rusanov_flux(name, wave_speed_default) result(this)
        character(len=*), optional, intent(in) :: name
        real(real64), optional, intent(in) :: wave_speed_default
        
        if (present(name)) then
            this%name = name
        else
            this%name = "Rusanov"
        end if
        
        if (present(wave_speed_default)) then
            this%wave_speed_default = wave_speed_default
        else
            this%wave_speed_default = 1.0_real64
        end if
    end function create_rusanov_flux
    
    subroutine rusanov_info(this)
        class(rusanov_flux), intent(in) :: this
        call flux_info(this)
        print *, "  Type: Rusanov flux"
        print *, "  Default wave speed: ", this%wave_speed_default
    end subroutine rusanov_info
    
end module rusanov_flux_module