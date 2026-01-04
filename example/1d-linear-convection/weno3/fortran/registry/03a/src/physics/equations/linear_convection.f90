! src/physics/equations/linear_convection.f90
module linear_convection_equation
    use precision_module, only: wp, ip
    use physics_interface, only: physics_equation
    implicit none
    private
    
    ! 具体方程类型 - 先声明
    type, extends(physics_equation) :: linear_convection_eq
        real(wp) :: wave_speed = 1.0_wp
    contains
        procedure :: flux => lc_flux
        procedure :: speed => lc_speed
    end type linear_convection_eq
    
    ! 公开接口
    public :: wp, ip
    public :: linear_convection_eq, create_linear_convection_eq
    
contains
    
    ! 构造函数
    function create_linear_convection_eq(wave_speed) result(eq)
        real(wp), intent(in), optional :: wave_speed
        type(linear_convection_eq) :: eq
        
        eq%name = "Linear Convection"
        if (present(wave_speed)) then
            eq%wave_speed = wave_speed
        else
            eq%wave_speed = 1.0_wp
        end if
    end function create_linear_convection_eq
    
    ! 方法实现
    pure function lc_flux(this, u) result(f)
        class(linear_convection_eq), intent(in) :: this
        real(wp), intent(in) :: u
        real(wp) :: f
        f = this%wave_speed * u
    end function lc_flux
    
    pure function lc_speed(this) result(a)
        class(linear_convection_eq), intent(in) :: this
        real(wp) :: a
        a = this%wave_speed
    end function lc_speed
    
end module linear_convection_equation