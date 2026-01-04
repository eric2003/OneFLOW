! src/numerics/reconstructor/base.f90
module reconstructor_base_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    public :: wp, reconstructor_base
    
    ! Base reconstructor type
    type, abstract :: reconstructor_base
        integer :: order = 1
        character(len=20) :: name = "Base"
    contains
        procedure(reconstruct_interface), deferred :: reconstruct
        procedure :: info => reconstructor_info
    end type reconstructor_base
    
    abstract interface
        subroutine reconstruct_interface(this, q, qL, qR)
            import :: reconstructor_base, wp
            class(reconstructor_base), intent(in) :: this
            real(wp), intent(in) :: q(:)
            real(wp), intent(out) :: qL(:)
            real(wp), intent(out) :: qR(:)
        end subroutine reconstruct_interface
    end interface
    
contains
    
    subroutine reconstructor_info(this)
        class(reconstructor_base), intent(in) :: this
        
        print *, "Reconstructor information:"
        print *, "  Name: ", trim(this%name)
        print *, "  Order: ", this%order
    end subroutine reconstructor_info
    
end module reconstructor_base_module