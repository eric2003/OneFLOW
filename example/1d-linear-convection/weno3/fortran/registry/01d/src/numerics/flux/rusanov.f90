! src/numerics/flux/rusanov.f90
module rusanov_flux_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use flux_base_module, only: flux_calculator_base, wp
    use registry_module, only: register_component_with_factory
    implicit none
    
    private
    public :: wp, rusanov_flux, create_rusanov
    
    type, extends(flux_calculator_base) :: rusanov_flux
        real(wp) :: wave_speed_default = 1.0_wp
    contains
        procedure :: compute => rusanov_compute
        procedure :: info => rusanov_info
    end type rusanov_flux
    
contains
    
    ! Factory function
    subroutine create_rusanov(instance)
        class(*), allocatable, intent(out) :: instance
        type(rusanov_flux), allocatable :: flux
        
        allocate(flux)
        flux%name = "Rusanov"
        flux%wave_speed_default = 1.0_wp
        
        call move_alloc(flux, instance)
    end subroutine create_rusanov
    
    ! Rusanov flux computation for linear convection
    subroutine rusanov_compute(this, qL, qR, flux, wave_speed)
        class(rusanov_flux), intent(in) :: this
        real(wp), intent(in) :: qL(:), qR(:)
        real(wp), intent(out) :: flux(:)
        real(wp), intent(in) :: wave_speed
        
        integer :: i, n
        real(wp) :: max_speed, fL, fR
        
        n = size(qL)
        
        if (size(qR) /= n .or. size(flux) /= n) then
            print *, "[ERROR] Array size mismatch in rusanov_compute"
            return
        end if
        
        ! Maximum wave speed
        max_speed = abs(wave_speed)
        
        ! Compute Rusanov flux at each interface
        do i = 1, n
            ! Linear convection flux: f(u) = wave_speed * u
            fL = wave_speed * qL(i)
            fR = wave_speed * qR(i)
            
            ! Rusanov flux formula
            flux(i) = 0.5_wp * (fL + fR) - 0.5_wp * max_speed * (qR(i) - qL(i))
        end do
    end subroutine rusanov_compute
    
    subroutine rusanov_info(this)
        class(rusanov_flux), intent(in) :: this
        call flux_info(this)
        print *, "  Type: Rusanov flux"
        print *, "  Default wave speed: ", this%wave_speed_default
    end subroutine rusanov_info
    
end module rusanov_flux_module