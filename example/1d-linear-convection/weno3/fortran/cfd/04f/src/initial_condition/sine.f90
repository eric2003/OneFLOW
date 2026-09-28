! ==================== sine.f90 ====================
! Sine wave initial condition for OneFLOW-CFD
!
! Purpose: Sine wave initial condition
! Author: OneFLOW-CFD Team
! Date: Created
! ===================================================

module sine_module
    use kinds, only: dp
    implicit none
    
    private
    public :: sine_initial_condition
    
    ! 常量定义
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    
contains
    
    ! Initial condition: sine wave
    function sine_initial_condition(x, L) result(u0)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: L  ! 域长度，可选参数
        real(dp) :: u0
        
        real(dp) :: domain_length
        
        ! 处理可选参数
        if (present(L)) then
            domain_length = L
        else
            domain_length = 2.0_dp  ! 默认值
        end if
        
        ! 计算正弦波
        u0 = sin(2.0_dp * PI * x / domain_length)
    end function sine_initial_condition
    
end module sine_module