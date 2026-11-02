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

    ! Add构造函数接口
    interface eno_reconstructor
        module procedure create_eno_reconstructor
    end interface

contains

    ! 构造函数
    type(eno_reconstructor) function create_eno_reconstructor(name, order, epsilon) result(this)
        character(len=*), optional, intent(in) :: name
        integer, optional, intent(in) :: order
        real(real64), optional, intent(in) :: epsilon

        ! Set默认值
        if (present(name)) then
            this%name = name
        else
            this%name = "ENO"
        end if

        if (present(order)) then
            this%order = order
        else
            this%order = 3
        end if

        if (present(epsilon)) then
            this%epsilon = epsilon
        else
            this%epsilon = 1.0e-6_real64
        end if
    end function create_eno_reconstructor

    subroutine eno_info(this)
        class(eno_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  Type: ENO reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine eno_info

end module eno_reconstructor_module