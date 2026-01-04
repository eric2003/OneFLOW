! src/numerics/flux/rusanov.f90
module rusanov_flux_module
    use, intrinsic :: iso_fortran_env, only: real64
    use flux_base_module, only: flux_calculator_base
    implicit none
    
    private
    public :: real64, rusanov_flux, create_rusanov_flux
    
    type, extends(flux_calculator_base) :: rusanov_flux
        real(real64) :: wave_speed_default = 1.0_real64
    contains
        procedure :: info => rusanov_info
    end type rusanov_flux
    
    ! 构造函数接口
    interface rusanov_flux
        module procedure create_rusanov_flux
    end interface
    
contains
    
    ! 构造函数
    type(rusanov_flux) function create_rusanov_flux() result(this)
        this%name = "Rusanov"
        this%wave_speed_default = 1.0_real64
    end function create_rusanov_flux
    
    subroutine rusanov_info(this)
        class(rusanov_flux), intent(in) :: this
        
        ! 调用父类的基础信息打印
        print *, "Flux calculator information:"
        call this%print_basic_info()
        
        ! 添加Rusanov特有信息
        print *, "  Type: Rusanov flux"
        print *, "  Default wave speed: ", this%wave_speed_default
    end subroutine rusanov_info
    
end module rusanov_flux_module