! src/numerics/flux/base.f90
module flux_base_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    public :: wp, flux_calculator_base
    
    ! Base flux calculator type
    type, abstract :: flux_calculator_base
        character(len=20) :: name = "Base"
    contains
        procedure(compute_interface), deferred :: compute
        procedure :: info => flux_info
    end type flux_calculator_base
    
    abstract interface
        subroutine compute_interface(this, qL, qR, flux, wave_speed)
            import :: flux_calculator_base, wp
            class(flux_calculator_base), intent(in) :: this
            real(wp), intent(in) :: qL(:), qR(:)
            real(wp), intent(out) :: flux(:)
            real(wp), intent(in) :: wave_speed
        end subroutine compute_interface
    end interface
    
contains
    
    subroutine flux_info(this)
        class(flux_calculator_base), intent(in) :: this
        print *, "Flux calculator information:"
        print *, "  Name: ", trim(this%name)
    end subroutine flux_info
    
end module flux_base_module