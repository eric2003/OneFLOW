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

    ! Add构造函数接口
    interface weno3_reconstructor
        module procedure create_weno3_reconstructor
    end interface

contains

    ! 构造函数
    type(weno3_reconstructor) function create_weno3_reconstructor(name, order, epsilon) result(this)
        character(len=*), optional, intent(in) :: name
        integer, optional, intent(in) :: order
        real(real64), optional, intent(in) :: epsilon

        if (present(name)) then
            this%name = name
        else
            this%name = "WENO3"
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
    end function create_weno3_reconstructor

    subroutine weno3_info(this)
        class(weno3_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  Type: WENO-3 reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine weno3_info

end module weno3_reconstructor_module