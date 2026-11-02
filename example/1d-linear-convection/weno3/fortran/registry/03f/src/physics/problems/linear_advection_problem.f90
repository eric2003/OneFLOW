! src/physics/problems/linear_advection_problem.f90
module linear_advection_problem
    use base_modules, only: wp, ip
    use initial_condition_module, only: create_initial_condition
    use boundary_condition_module, only: create_boundary_condition
    implicit none
    private
    
    type, public :: linear_advection_problem
        character(len=20) :: ic_type = "step"
        character(len=20) :: boundary_type = "periodic"
        real(wp) :: left_value = 1.0_wp
        real(wp) :: right_value = 2.0_wp
        real(wp) :: domain_length = 2.0_wp
        real(wp) :: wave_speed = 1.0_wp
    contains
        procedure :: create_ic => prob_create_ic
        procedure :: create_bc => prob_create_bc
        procedure :: exact_solution => prob_exact_solution
    end type
    
contains

    function prob_create_ic(this, config) result(ic)
        class(linear_advection_problem), intent(in) :: this
        class(*), intent(in) :: config
        class(*), allocatable :: ic
        
        ! 通过初始条件工厂创建
        call create_initial_condition(this%ic_type, config, ic)
    end function
    
    function prob_create_bc(this, cfd) result(bc)
        class(linear_advection_problem), intent(in) :: this
        class(*), intent(in) :: cfd
        class(*), allocatable :: bc
        
        ! 通过边界条件工厂创建
        call create_boundary_condition(this%boundary_type, cfd, bc)
    end function
    
    function prob_exact_solution(this, x, t) result(u)
        class(linear_advection_problem), intent(in) :: this
        real(wp), intent(in) :: x(:), t
        real(wp), allocatable :: u(:)
        
        integer :: i, n
        real(wp) :: x_shifted, L
        
        n = size(x)
        L = this%domain_length
        allocate(u(n))
        
        ! 周期性平移
        do i = 1, n
            x_shifted = x(i) - this%wave_speed * t
            x_shifted = modulo(x_shifted, L)
            if (x_shifted < 0.0_wp) x_shifted = x_shifted + L
            
            ! 初始条件逻辑（简化）
            if (this%ic_type == "step" .and. &
                x_shifted >= 0.5_wp .and. x_shifted <= 1.0_wp) then
                u(i) = 2.0_wp
            else
                u(i) = 1.0_wp
            end if
        end do
    end function
    
end module