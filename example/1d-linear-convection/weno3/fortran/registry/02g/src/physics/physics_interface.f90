! src/physics/physics_interface.f90
module physics_interface
    use precision_module, only: wp, ip
    implicit none
    private
    
    ! 定义抽象基类型
    type, abstract :: physics_equation
        character(len=:), allocatable :: name
    contains
        procedure(eq_flux_abs), deferred :: flux
        procedure(eq_speed_abs), deferred :: speed
    end type physics_equation
    
    type, abstract :: physics_problem
        character(len=:), allocatable :: name
    contains
        procedure(prob_ic_abs), deferred :: initial_condition
        procedure(prob_bc_abs), deferred :: boundary_condition
        procedure(prob_exact_abs), deferred :: exact_solution
    end type physics_problem
    
    ! 抽象接口定义
    abstract interface
        pure function eq_flux_abs(this, u) result(f)
            import :: physics_equation, wp
            class(physics_equation), intent(in) :: this
            real(wp), intent(in) :: u
            real(wp) :: f
        end function eq_flux_abs
        
        pure function eq_speed_abs(this) result(a)
            import :: physics_equation, wp
            class(physics_equation), intent(in) :: this
            real(wp) :: a
        end function eq_speed_abs
        
        subroutine prob_ic_abs(this, x, u)
            import :: physics_problem, wp
            class(physics_problem), intent(in) :: this
            real(wp), intent(in) :: x(:)
            real(wp), intent(out) :: u(:)
        end subroutine prob_ic_abs
        
        subroutine prob_bc_abs(this, u, t)
            import :: physics_problem, wp
            class(physics_problem), intent(in) :: this
            real(wp), intent(inout) :: u(:)
            real(wp), intent(in), optional :: t
        end subroutine prob_bc_abs
        
        function prob_exact_abs(this, x, t) result(u)
            import :: physics_problem, wp
            class(physics_problem), intent(in) :: this
            real(wp), intent(in) :: x(:), t
            real(wp), dimension(size(x)) :: u
        end function prob_exact_abs
    end interface
    
    ! 公开接口
    public :: physics_equation, physics_problem
    public :: wp, ip
    
end module physics_interface