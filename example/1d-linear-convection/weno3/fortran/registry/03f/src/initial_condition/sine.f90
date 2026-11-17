! src/initial_condition/sine.f90
module sine_ic_module
    use base_modules, only: wp, ip
    use ic_base_module, only: initial_condition
    use solution_module, only: solution_type
    use domain_module, only: domain_type
    implicit none
    private
    
    type, extends(initial_condition), public :: sine_wave_ic
    contains
        procedure :: evaluate_at => sine_evaluate_at
        procedure :: apply => sine_apply
    end type
    
    interface sine_wave_ic
        module procedure create_sine_wave_ic
    end interface
    
contains

    type(sine_wave_ic) function create_sine_wave_ic() result(this)
        this%name = "sin"
    end function
    
    function sine_evaluate_at(this, x) result(u)
        class(sine_wave_ic), intent(in) :: this
        real(wp), intent(in) :: x(:)
        real(wp) :: u(size(x))
        
        integer :: i
        real(wp) :: L
        
        ! 假设域长度，可以根据需要调整
        L = 2.0_wp  ! 默认域长度
        
        do i = 1, size(x)
            u(i) = sin(2.0_wp * 3.141592653589793_wp * x(i) / L)
        end do
    end function
    
    subroutine sine_apply(this, solution)
        class(sine_wave_ic), intent(in) :: this
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
    
end module sine_ic_module