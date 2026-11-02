! ==================== rusanov_flux.f90 ====================
! Rusanov flux computation for OneFLOW-CFD
!
! Purpose: Compute numerical flux using Rusanov scheme
! Reference: V.V. Rusanov, "Calculation of Interaction of Non-steady
!            Shock Waves with Obstacles," J. Comp. Math. Phys. USSR, 1961
! Author: OneFLOW-CFD Team
! Date: Created
! =========================================================

module rusanov_flux_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    implicit none
    
    private
    public :: rusanov_flux
    
contains
    
    ! Rusanov flux implementation
    subroutine rusanov_flux(q_face_left, q_face_right, flux, config, mesh)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdConfigType), intent(in) :: config
        type(MeshType), intent(in) :: mesh
        
        integer :: i
        real(dp) :: u_L, u_R, F_L, F_R, c_L, c_R, Smax
        
        c_L = config%wave_speed
        c_R = config%wave_speed
        
        do i = 1, mesh%nnodes
            u_L = q_face_left(i)
            u_R = q_face_right(i)
            F_L = c_L * u_L
            F_R = c_R * u_R
            Smax = max(abs(c_L), abs(c_R))
            flux(i) = 0.5_dp * (F_L + F_R) - 0.5_dp * Smax * (u_R - u_L)
        end do
    end subroutine rusanov_flux
    
end module rusanov_flux_module