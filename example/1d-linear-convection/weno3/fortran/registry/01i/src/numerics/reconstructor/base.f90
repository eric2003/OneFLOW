! src/numerics/reconstructor/base.f90
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
        procedure :: reconstruct => reconstruct_default  ! 添加这个方法
    end type reconstructor_base

contains

    subroutine reconstructor_info(this)
        class(reconstructor_base), intent(in) :: this
        print *, "Reconstructor information:"
        print *, "  Name: ", trim(this%name)
        print *, "  Order: ", this%order
    end subroutine reconstructor_info

    ! 默认的reconstructionmethod
    subroutine reconstruct_default(this, q, qL, qR)
        class(reconstructor_base), intent(in) :: this
        real(real64), intent(in) :: q(0:)  ! 包含ghost cells
        real(real64), intent(out) :: qL(:), qR(:)
        integer :: i, n

        n = size(qL)
        do i = 1, n
            qL(i) = q(i)        ! 简单的一阶重构
            qR(i) = q(i+1)
        end do
    end subroutine reconstruct_default
end module reconstructor_base_module