! src/numerics/reconstructor/eno.f90
module eno_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: real64
    use reconstructor_base_module, only: reconstructor_base
    implicit none
    
    private
	public :: real64, eno_reconstructor, create_eno_reconstructor  ! ← 添加这个
    
    type, extends(reconstructor_base) :: eno_reconstructor
        real(real64) :: epsilon = 1.0e-6_real64
    contains
        procedure :: info => eno_info
    end type eno_reconstructor
    
    ! 构造函数接口
    interface eno_reconstructor
        module procedure create_eno_reconstructor
    end interface
    
contains
    
    ! 构造函数
    type(eno_reconstructor) function create_eno_reconstructor() result(this)
        this%name = "ENO"
        this%order = 3
        this%epsilon = 1.0e-6_real64
    end function create_eno_reconstructor
    
    subroutine eno_info(this)
        class(eno_reconstructor), intent(in) :: this
        
        ! 调用父类的基础信息打印
        print *, "Reconstructor information:"
        call this%print_basic_info()
        
        ! 添加ENO特有信息
        print *, "  Type: ENO reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine eno_info
    
end module eno_reconstructor_module