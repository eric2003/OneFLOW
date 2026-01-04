! src/initial_condition/gaussian.f90
module gaussian_ic_module
    use base_modules, only: wp, ip
    use ic_base_module, only: initial_condition
    use solution_module, only: solution_type
    use domain_module, only: domain_type
    implicit none
    private
    
    type, extends(initial_condition), public :: gaussian_pulse_ic
    contains
        procedure :: evaluate_at => gaussian_evaluate_at
        procedure :: apply => gaussian_apply
    end type
    
    interface gaussian_pulse_ic
        module procedure create_gaussian_pulse_ic
    end interface
    
contains

    type(gaussian_pulse_ic) function create_gaussian_pulse_ic() result(this)
        this%name = "gaussian"
    end function
    
    function gaussian_evaluate_at(this, x) result(u)
        class(gaussian_pulse_ic), intent(in) :: this
        real(wp), intent(in) :: x(:)
        real(wp) :: u(size(x))
        
        integer :: i
        real(wp) :: center, width
        
        center = 0.5_wp  ! 脉冲中心
        width = 0.1_wp   ! 脉冲宽度
        
        do i = 1, size(x)
            u(i) = exp(-((x(i) - center) / width)**2)
        end do
    end function
    
    subroutine gaussian_apply(this, solution)
        class(gaussian_pulse_ic), intent(in) :: this
        type(solution_type), intent(inout) :: solution
        
        integer :: i, idx
        real(wp), allocatable :: u0(:)
        type(domain_type), pointer :: domain
        
        domain => solution%domain
        
        ! 评估初始条件
        u0 = this%evaluate_at(domain%mesh%xcc)
        
        ! 应用到物理区域
        do i = domain%ist, domain%ied - 1
            idx = i - domain%ist + 1
            solution%u(i) = u0(idx)
        end do
    end subroutine
    
end module gaussian_ic_module