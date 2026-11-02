! ==================== step.f90 =====================
! Step function initial condition for OneFLOW-CFD
!
! Purpose: Step function (Riemann problem) initial condition
! Author: OneFLOW-CFD Team
! Date: Created
! ===================================================

module step_module
    use kinds, only: dp
    implicit none
    
    private
    public :: step_initial_condition
    
contains
    
    ! Initial condition: step function
    function step_initial_condition(x) result(u0)
        real(dp), intent(in) :: x
        real(dp) :: u0
        
        if (0.5_dp <= x .and. x <= 1.0_dp) then
            u0 = 2.0_dp
        else
            u0 = 1.0_dp
        end if
    end function step_initial_condition
    
end module step_module