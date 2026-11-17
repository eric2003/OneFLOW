! ==================== initial_condition.f90 ==========================
! Initial conditions and analytical solutions for OneFLOW-CFD
!
! Purpose: Main module for initial conditions and analytical solutions
! Author: OneFLOW-CFD Team
! Date: Created
! ======================================================================

module initial_condition_module
    use kinds, only: dp
    use config_module, only: CfdConfigType

    use step_module, only: step_initial_condition
    use sine_module, only: sine_initial_condition	
    
    implicit none
    
    private
    public :: initial_condition, analytical_solution, init_ic
	
	type(CfdConfigType) :: config
	
    integer, parameter :: IC_STEP = 1
    integer, parameter :: IC_SINE = 2
    
    integer, save :: current_ic_type = IC_STEP
    real(dp), save :: domain_length = 2.0_dp  ! 用于正弦波的域长度
    
contains
    subroutine init_ic(ic_type, L)
        character(len=10), intent(in) :: ic_type
        real(dp), intent(in), optional :: L
        
        if (present(L)) then
            domain_length = L
        end if
        
        select case(ic_type)
        case("step")
            current_ic_type = IC_STEP
        case("sine")
            current_ic_type = IC_SINE
        case default
            current_ic_type = IC_STEP
        end select
    end subroutine init_ic
	
    function initial_condition(x) result(u0)
        real(dp), intent(in) :: x
        real(dp) :: u0
        
        select case(current_ic_type)
        case(IC_STEP)
            u0 = step_initial_condition(x)
        case(IC_SINE)
            u0 = sine_initial_condition(x, domain_length)
        case default
            u0 = step_initial_condition(x)
        end select
    end function initial_condition
    
    ! Analytical solution with periodic BC (保持原有功能)
    function analytical_solution(x, t, a, L) result(u)
        real(dp), intent(in) :: x, t, a, L
        real(dp) :: u, x_shifted
        
        x_shifted = mod(x - a * t + L, L)
        u = initial_condition(x_shifted)
    end function analytical_solution
    
end module initial_condition_module