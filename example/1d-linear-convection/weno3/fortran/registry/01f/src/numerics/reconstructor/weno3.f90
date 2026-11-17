! src/numerics/reconstructor/weno3.f90
module weno3_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: real64
    use reconstructor_base_module, only: reconstructor_base
    implicit none
    
    private
    public :: real64, weno3_reconstructor
    
    type, extends(reconstructor_base) :: weno3_reconstructor
        real(real64) :: epsilon = 1.0e-6_real64
    contains
        procedure :: info => weno3_info
    end type weno3_reconstructor
    
contains
    
    subroutine weno3_info(this)
        class(weno3_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  Type: WENO-3 reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine weno3_info
    
end module weno3_reconstructor_module