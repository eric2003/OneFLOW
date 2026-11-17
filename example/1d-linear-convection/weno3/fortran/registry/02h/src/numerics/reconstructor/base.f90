!src/numerics/reconstructor/base.f90
module reconstructor_base_module
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    
    private
    public :: real64, reconstructor_base
    
    ! Base reconstructor type
    type :: reconstructor_base
        integer :: order = 1
        character(len=20) :: name = "Base"
    contains
        procedure :: info => reconstructor_info
        procedure :: print_basic_info => reconstructor_print_basic  ! 添加一个辅助方法
    end type reconstructor_base
    
contains
    
    subroutine reconstructor_info(this)
        class(reconstructor_base), intent(in) :: this
        print *, "Reconstructor information:"
        print *, "  Name: ", trim(this%name)
        print *, "  Order: ", this%order
    end subroutine reconstructor_info
    
    subroutine reconstructor_print_basic(this)
        class(reconstructor_base), intent(in) :: this
        print *, "  Name: ", trim(this%name)
        print *, "  Order: ", this%order
    end subroutine reconstructor_print_basic
    
end module reconstructor_base_module