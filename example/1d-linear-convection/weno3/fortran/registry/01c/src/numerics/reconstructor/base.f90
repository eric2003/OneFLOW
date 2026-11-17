! src/numerics/reconstructor/base.f90
module reconstructor_base_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    public :: wp
    
    ! Base reconstructor type (simplified)
    type :: reconstructor_base
        integer :: order = 1
        character(len=20) :: name = "Base"
    contains
        procedure :: info => reconstructor_info
    end type reconstructor_base
    
contains
    
    subroutine reconstructor_info(this)
        class(reconstructor_base), intent(in) :: this
        
        print *, "Reconstructor information:"
        print *, "  Name: ", trim(this%name)
        print *, "  Order: ", this%order
    end subroutine reconstructor_info
    
end module reconstructor_base_module