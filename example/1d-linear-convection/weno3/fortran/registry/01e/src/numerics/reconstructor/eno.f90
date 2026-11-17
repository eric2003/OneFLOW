! src/numerics/reconstructor/eno.f90
module eno_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: real64
    use reconstructor_base_module, only: reconstructor_base
    implicit none
    
    private
    public :: real64, eno_reconstructor
    
    type, extends(reconstructor_base) :: eno_reconstructor
        real(real64) :: epsilon = 1.0e-6_real64
    contains
        procedure :: info => eno_info
    end type eno_reconstructor
    
contains
    
    subroutine eno_info(this)
        class(eno_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  Type: ENO reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine eno_info
    
end module eno_reconstructor_module