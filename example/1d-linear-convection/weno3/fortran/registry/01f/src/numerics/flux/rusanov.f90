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
    
contains
    
    subroutine rusanov_info(this)
        class(rusanov_flux), intent(in) :: this
        call flux_info(this)
        print *, "  Type: Rusanov flux"
        print *, "  Default wave speed: ", this%wave_speed_default
    end subroutine rusanov_info
    
end module rusanov_flux_module