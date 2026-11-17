! src/physics/equations/linear_convection.f90
module linear_convection_equation
    use precision_module, only: wp, ip
    use physics_interface, only: physics_equation
    implicit none
    private
    
    public :: linear_convection_eq, create_linear_convection_eq
    public :: linear_convection_flux, linear_convection_speed
    
    ! 具体方程类型
    type, extends(physics_equation) :: linear_convection_eq
        real(wp) :: wave_speed = 1.0_wp
    contains
        procedure :: flux => lc_flux
        procedure :: speed => lc_speed
    end type linear_convection_eq
    
    ! 独立函数（供外部调用）
    interface linear_convection_flux
        module procedure :: lc_flux_func
    end interface
    
    interface linear_convection_speed
        module procedure :: lc_speed_func
    end interface
    
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
    
    ! 独立函数
    pure function lc_flux_func(u, wave_speed) result(f)
        real(wp), intent(in) :: u
        real(wp), intent(in), optional :: wave_speed
        real(wp) :: f
        real(wp) :: a
        
        if (present(wave_speed)) then
            a = wave_speed
        else
            a = 1.0_wp
        end if
        f = a * u
    end function lc_flux_func
    
    pure function lc_speed_func(wave_speed) result(a)
        real(wp), intent(in), optional :: wave_speed
        real(wp) :: a
        
        if (present(wave_speed)) then
            a = wave_speed
        else
            a = 1.0_wp
        end if
    end function lc_speed_func
    
end module linear_convection_equation