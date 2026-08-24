! ==================== initial_condition.f90 ====================
! Initial conditions and analytical solutions for OneFLOW-CFD

module initial_condition_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    private
    
    public :: initial_condition, analytical_solution
    
contains
    
    ! ===================================================================
    ! Initial Conditions
    ! ===================================================================
    
    ! Initial condition: step function
    function initial_condition(x) result(u0)
        real(dp), intent(in) :: x
        real(dp) :: u0
        
        if (0.5_dp <= x .and. x <= 1.0_dp) then
            u0 = 2.0_dp
        else
            u0 = 1.0_dp
        end if
    end function initial_condition
    
    ! Analytical solution with periodic BC
    function analytical_solution(x, t, a, L) result(u)
        real(dp), intent(in) :: x, t, a, L
        real(dp) :: u, x_shifted
        
        x_shifted = mod(x - a * t + L, L)
        u = initial_condition(x_shifted)
    end function analytical_solution
    
end module initial_condition_module