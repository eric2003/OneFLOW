! src/physics/equations/linear_advection.f90
module linear_advection_equation
    use base_modules, only: wp, ip
    implicit none
    private
    
    type, public :: linear_advection_eq
        real(wp) :: wave_speed = 1.0_wp
    contains
        procedure :: flux => eq_flux
        procedure :: max_wave_speed => eq_max_wave_speed
        procedure :: num_eqs => eq_num_eqs
        procedure :: print_info => eq_print_info
    end type
    
contains

    pure function eq_flux(this, u) result(f)
        class(linear_advection_eq), intent(in) :: this
        real(wp), intent(in) :: u
        real(wp) :: f
        f = this%wave_speed * u
    end function
    
    pure function eq_max_wave_speed(this, u) result(smax)
        class(linear_advection_eq), intent(in) :: this
        real(wp), intent(in) :: u
        real(wp) :: smax
        smax = abs(this%wave_speed)
    end function
    
    integer function eq_num_eqs(this)
        class(linear_advection_eq), intent(in) :: this
        eq_num_eqs = 1
    end function
    
    subroutine eq_print_info(this)
        class(linear_advection_eq), intent(in) :: this
        print *, "Linear Advection Equation:"
        print *, "  Wave speed: ", this%wave_speed
    end subroutine
    
end module