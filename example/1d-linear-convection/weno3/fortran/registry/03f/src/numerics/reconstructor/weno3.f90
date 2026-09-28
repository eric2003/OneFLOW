! src/numerics/reconstructor/weno3.f90
module weno3_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: real64
    use reconstructor_base_module, only: reconstructor_base
    implicit none
    private
    
    type, extends(reconstructor_base), public :: weno3_reconstructor
        real(real64) :: epsilon = 1.0e-6_real64
    contains
        procedure :: info => weno3_info
    end type weno3_reconstructor
    
    ! 构造函数接口 - 使用类型名本身作为构造函数
    interface weno3_reconstructor
        module procedure create_weno3_reconstructor
    end interface
    
contains
    
    ! 构造函数
    type(weno3_reconstructor) function create_weno3_reconstructor() result(this)
        this%name = "WENO3"
        this%order = 3
        this%epsilon = 1.0e-6_real64
    end function create_weno3_reconstructor
    
    subroutine weno3_info(this)
        class(weno3_reconstructor), intent(in) :: this
        
        ! 调用父类的基础信息打印
        print *, "Reconstructor information:"
        call this%print_basic_info()
        
        ! 添加WENO3特有信息
        print *, "  Type: WENO-3 reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine weno3_info
    
end module weno3_reconstructor_module