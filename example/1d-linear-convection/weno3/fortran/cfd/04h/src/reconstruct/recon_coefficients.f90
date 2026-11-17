! ==================== recon_coefficients.f90 ====================
! Reconstruction coefficients for OneFLOW-CFD
!
! Purpose: Initialize coefficients for ENO/WENO reconstruction schemes
! Author: OneFLOW-CFD Team
! Date: Created
! ===============================================================

module recon_coefficients_module
    use kinds, only: dp
    implicit none
    
    private
    public :: init_coef
    
contains
    
    ! Initialize ENO/WENO coefficients
    subroutine init_coef(spatial_order, coef)
        integer, intent(in) :: spatial_order
        real(dp), intent(out) :: coef(:,:)
        
        coef = 0.0_dp
        
        select case(spatial_order)
        case(1)
            coef(1,1) = 1.0_dp
            coef(2,1) = 1.0_dp
            
        case(2)
            coef(1,1:2) = [ 3.0_dp/2.0_dp, -1.0_dp/2.0_dp ]
            coef(2,1:2) = [ 1.0_dp/2.0_dp,  1.0_dp/2.0_dp ]
            coef(3,1:2) = [ -1.0_dp/2.0_dp, 3.0_dp/2.0_dp ]
            
        case(3)
            coef(1,1:3) = [ 11.0_dp/6.0_dp, -7.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
            coef(2,1:3) = [  1.0_dp/3.0_dp,  5.0_dp/6.0_dp, -1.0_dp/6.0_dp ]
            coef(3,1:3) = [ -1.0_dp/6.0_dp,  5.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
            coef(4,1:3) = [  1.0_dp/3.0_dp, -7.0_dp/6.0_dp, 11.0_dp/6.0_dp ]
            
        case default
            error stop "Unsupported spatial order"
        end select
    end subroutine init_coef
    
end module recon_coefficients_module