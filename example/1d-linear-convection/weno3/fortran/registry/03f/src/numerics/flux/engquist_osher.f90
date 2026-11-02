! src/numerics/flux/engquist_osher.f90
module engquist_osher_flux
    use base_modules, only: wp, ip
    use flux_base_module, only: flux_calculator_base
    implicit none
    private
    
    type, extends(flux_calculator_base), public :: engquist_osher_flux_t
        real(wp) :: wave_speed = 1.0_wp
    contains
        procedure :: compute => eo_compute
    end type
    
contains

    subroutine eo_compute(this, qL, qR, flux, equation)
        class(engquist_osher_flux_t), intent(inout) :: this
        real(wp), intent(in) :: qL(:), qR(:)
        real(wp), intent(out) :: flux(:)
        class(*), intent(in) :: equation
        
        integer :: i, n
        real(wp) :: cp, cm, uL, uR
        
        select type(eq => equation)
        type is (linear_advection_eq)
            cp = 0.5_wp * (eq%wave_speed + abs(eq%wave_speed))
            cm = 0.5_wp * (eq%wave_speed - abs(eq%wave_speed))
            
            n = size(qL)
            do i = 1, n
                uL = qL(i)
                uR = qR(i)
                flux(i) = cp * uL + cm * uR
            end do
        class default
            error stop "Unsupported equation type for Engquist-Osher flux"
        end select
    end subroutine
    
end module