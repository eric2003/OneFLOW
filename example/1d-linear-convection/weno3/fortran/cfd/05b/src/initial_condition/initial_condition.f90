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
    use complex_wave_precise_module, only: complex_wave_precise, init_precise_params
    implicit none
    
    private
    public :: initial_condition, analytical_solution, init_ic
	
	type(CfdConfigType) :: config
	
    integer, parameter :: IC_STEP = 1
    integer, parameter :: IC_SINE = 2
	integer, parameter :: IC_COMPLEX = 8  ! 新增复杂波形
    
    integer, save :: current_ic_type = IC_STEP
    real(dp), save :: domain_length = 2.0_dp  ! 用于正弦波的域长度
    
contains
    subroutine init_ic(ic_type, L)
        character(len=10), intent(in) :: ic_type
        real(dp), intent(in), optional :: L
		real(dp) :: beta, z, alpha, a, delta
        
        if (present(L)) then
            domain_length = L
        end if
        
        select case(ic_type)
        case("step")
            current_ic_type = IC_STEP
        case("sine")
            current_ic_type = IC_SINE
        case("complex")
            current_ic_type = IC_COMPLEX
            call init_precise_params()
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
        case(IC_COMPLEX)
            u0 = complex_wave_precise(x)
        case default
            u0 = step_initial_condition(x)
        end select
    end function initial_condition
    
    ! Analytical solution with periodic BC (保持原有功能)
    function analytical_solution(x, t, wave_speed, xmin, xmax) result(u)
        real(dp), intent(in) :: x, t, wave_speed, xmin, xmax
        real(dp) :: u, x_shifted
        
        ! 使用周期性边界平移坐标
        x_shifted = shift_with_periodic_bc(x, t * wave_speed, xmin, xmax)        
        
        u = initial_condition(x_shifted)
    end function analytical_solution
    
    ! ===================================================================
    ! 周期性边界下的坐标平移（辅助函数）
    ! ===================================================================
    pure function shift_with_periodic_bc(x, shift, xmin, xmax) result(x_shifted)
        real(dp), intent(in) :: x, shift, xmin, xmax
        real(dp) :: x_shifted, L
        
        L = xmax - xmin
        
        if (L <= 0.0_dp) then
            x_shifted = x  ! 无效域，返回原值
            return
        end if
        
        ! 应用平移
        x_shifted = x - shift
        
        ! 使用 Fortran 的 modulo 函数处理周期性边界
        ! modulo(a, p) 返回 a - floor(a/p) * p，适用于负数
        x_shifted = xmin + modulo(x_shifted - xmin, L)
        
        ! 处理边界情况（由于浮点误差）
        if (x_shifted >= xmax - epsilon(1.0_dp)) then
            x_shifted = xmin
        end if
    end function shift_with_periodic_bc    
    
    
end module initial_condition_module