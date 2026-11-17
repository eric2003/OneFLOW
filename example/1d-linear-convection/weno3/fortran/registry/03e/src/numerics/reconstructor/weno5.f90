! src/numerics/reconstructor/weno5.f90（新增）
module weno5_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: real64
    use reconstructor_base_module, only: reconstructor_base
    implicit none
    
    private
    public :: real64, weno5_reconstructor, create_weno5_reconstructor
    
    type, extends(reconstructor_base) :: weno5_reconstructor
        real(real64) :: epsilon = 1.0e-6_real64
    contains
        procedure :: info => weno5_info
    end type weno5_reconstructor
    
    ! 构造函数接口
    interface weno5_reconstructor
        module procedure create_weno5_reconstructor
    end interface
    
contains
    
    ! 构造函数
    type(weno5_reconstructor) function create_weno5_reconstructor() result(this)
        this%name = "WENO5"
        this%order = 5
        this%epsilon = 1.0e-6_real64
    end function create_weno5_reconstructor
    
    subroutine weno5_info(this)
        class(weno5_reconstructor), intent(in) :: this
        
        ! 调用父类的基础信息打印
        print *, "Reconstructor information:"
        call this%print_basic_info()
        
        ! 添加WENO5特有信息
        print *, "  Type: WENO-5 reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine weno5_info
    
end module weno5_reconstructor_module