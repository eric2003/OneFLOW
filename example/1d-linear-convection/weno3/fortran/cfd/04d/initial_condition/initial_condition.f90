! ==================== initial_condition.f90 ==========================
! Initial conditions and analytical solutions for OneFLOW-CFD
!
! Purpose: Main module for initial conditions and analytical solutions
! Author: OneFLOW-CFD Team
! Date: Created
! ======================================================================

module initial_condition_module
    use kinds, only: dp
    use step_module, only: step_initial_condition
    
    implicit none
    
    private
    public :: initial_condition, analytical_solution
    
    ! 保持原有函数名兼容性
    interface initial_condition
        module procedure step_initial_condition
    end interface
    
contains
    
    ! Analytical solution with periodic BC (保持原有功能)
    function analytical_solution(x, t, a, L) result(u)
        real(dp), intent(in) :: x, t, a, L
        real(dp) :: u, x_shifted
        
        x_shifted = mod(x - a * t + L, L)
        u = initial_condition(x_shifted)
    end function analytical_solution
    
end module initial_condition_module