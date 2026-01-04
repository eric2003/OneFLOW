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
    
    ! 必须添加构造函数接口
    interface eno_reconstructor
        module procedure create_eno_reconstructor
    end interface
    
contains
    
    ! 构造函数实现
    type(eno_reconstructor) function create_eno_reconstructor() result(this)
        this%name = "ENO"
        this%order = 3
        this%epsilon = 1.0e-6_real64
    end function create_eno_reconstructor
    
    subroutine eno_info(this)
        class(eno_reconstructor), intent(in) :: this
        ! 必须调用父类的info方法
        call reconstructor_info(this)  ! 这会调用base模块中的reconstructor_info
        print *, "  Type: ENO reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine eno_info
    
end module eno_reconstructor_module