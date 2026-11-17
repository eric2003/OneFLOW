! src/physics/problems/linear_convection_problem.f90
module linear_convection_problem
    use precision_module, only: wp, ip
    use physics_interface, only: physics_problem
    implicit none
    private
    
    public :: linear_convection_prob, create_linear_convection_prob
    
    ! 具体问题类型
    type, extends(physics_problem) :: linear_convection_prob
        real(wp) :: wave_speed = 1.0_wp
        real(wp) :: domain_length = 2.0_wp
        character(len=20) :: ic_type = "step"
        character(len=20) :: boundary_type = "periodic"
    contains
        procedure :: initial_condition => lc_initial_condition
        procedure :: boundary_condition => lc_boundary_condition
        procedure :: exact_solution => lc_exact_solution
    end type linear_convection_prob
    
contains
    
    ! 构造函数
    function create_linear_convection_prob(wave_speed, domain_length, &
                                          ic_type, boundary_type) result(prob)
        real(wp), intent(in), optional :: wave_speed, domain_length
        character(len=*), intent(in), optional :: ic_type, boundary_type
        type(linear_convection_prob) :: prob
        
        prob%name = "Linear Convection Problem"
        
        if (present(wave_speed)) prob%wave_speed = wave_speed
        if (present(domain_length)) prob%domain_length = domain_length
        if (present(ic_type)) prob%ic_type = ic_type
        if (present(boundary_type)) prob%boundary_type = boundary_type
    end function create_linear_convection_prob
    
    ! 初始条件
    subroutine lc_initial_condition(this, x, u)
        class(linear_convection_prob), intent(in) :: this
        real(wp), intent(in) :: x(:)
        real(wp), intent(out) :: u(:)
        
        integer :: i
        
        select case (trim(this%ic_type))
        case ("step")
            do i = 1, size(x)
                if (x(i) >= 0.5_wp .and. x(i) <= 1.0_wp) then
                    u(i) = 2.0_wp
                else
                    u(i) = 1.0_wp
                end if
            end do
            
        case ("sin", "sine")
            do i = 1, size(x)
                u(i) = sin(2.0_wp * 3.141592653589793_wp * x(i) / this%domain_length)
            end do
            
        case ("gaussian")
            do i = 1, size(x)
                u(i) = exp(-((x(i) - 0.5_wp) / 0.1_wp)**2)
            end do
            
        case default
            ! 默认阶跃函数
            do i = 1, size(x)
                if (x(i) >= 0.5_wp .and. x(i) <= 1.0_wp) then
                    u(i) = 2.0_wp
                else
                    u(i) = 1.0_wp
                end if
            end do
        end select
    end subroutine lc_initial_condition
    
    ! 边界条件（虚拟实现，实际在boundary模块）
    subroutine lc_boundary_condition(this, u, t)
        class(linear_convection_prob), intent(in) :: this
        real(wp), intent(inout) :: u(:)
        real(wp), intent(in), optional :: t
        
        ! 边界条件将在独立模块实现
        print *, "[PROBLEM] Boundary condition placeholder"
        if (present(t)) then
            print *, "  Time = ", t
        end if
    end subroutine lc_boundary_condition
    
    ! 精确解（周期性平移）
    function lc_exact_solution(this, x, t) result(u)
        class(linear_convection_prob), intent(in) :: this
        real(wp), intent(in) :: x(:), t
        real(wp), dimension(size(x)) :: u
        real(wp), dimension(size(x)) :: x_shifted
        integer :: i
        
        ! 周期性平移
        do i = 1, size(x)
            x_shifted(i) = x(i) - this%wave_speed * t
            ! 确保在 [0, domain_length) 范围内
            do while (x_shifted(i) < 0.0_wp)
                x_shifted(i) = x_shifted(i) + this%domain_length
            end do
            do while (x_shifted(i) >= this%domain_length)
                x_shifted(i) = x_shifted(i) - this%domain_length
            end do
        end do
        
        ! 重用初始条件函数
        call this%initial_condition(x_shifted, u)
    end function lc_exact_solution
    
end module linear_convection_problem