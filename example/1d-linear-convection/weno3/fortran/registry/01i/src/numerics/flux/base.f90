! src/numerics/flux/base.f90
module flux_base_module
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    
    private
    public :: real64, flux_calculator_base
    
    type :: flux_calculator_base
        character(len=20) :: name = "Base"
    contains
        procedure :: info => flux_info
    end type flux_calculator_base
    
contains
    
    subroutine flux_info(this)
        class(flux_calculator_base), intent(in) :: this
        print *, "Flux calculator information:"
        print *, "  Name: ", trim(this%name)
    end subroutine flux_info
    
end module flux_base_module