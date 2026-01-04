! src/initial_condition/step.f90
module step_ic_module
    use base_modules, only: wp, ip
    use ic_base_module, only: initial_condition
    use solution_module, only: solution_type
    use domain_module, only: domain_type
    implicit none
    private
    
    type, extends(initial_condition), public :: step_function_ic
    contains
        procedure :: evaluate_at => step_evaluate_at
        procedure :: apply => step_apply
    end type
    
    interface step_function_ic
        module procedure create_step_function_ic
    end interface
    
contains

    type(step_function_ic) function create_step_function_ic() result(this)
        this%name = "step"
    end function
    
    function step_evaluate_at(this, x) result(u)
        class(step_function_ic), intent(in) :: this
        real(wp), intent(in) :: x(:)
        real(wp) :: u(size(x))
        
        integer :: i
        
        do i = 1, size(x)
            if (x(i) >= 0.5_wp .and. x(i) <= 1.0_wp) then
                u(i) = 2.0_wp
            else
                u(i) = 1.0_wp
            end if
        end do
    end function
    
    subroutine step_apply(this, solution)
        class(step_function_ic), intent(in) :: this
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
    
end module step_ic_module