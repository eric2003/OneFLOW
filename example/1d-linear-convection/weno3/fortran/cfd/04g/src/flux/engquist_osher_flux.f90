! ==================== engquist_osher_flux.f90 ====================
! Engquist-Osher flux computation for OneFLOW-CFD
!
! Purpose: Compute numerical flux using Engquist-Osher scheme
! Reference: B. Engquist and S. Osher, "One-sided difference
!            approximations for nonlinear conservation laws,"
!            Math. Comp., 1981
! Author: OneFLOW-CFD Team
! Date: Created
! ================================================================

module engquist_osher_flux_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    implicit none
    
    private
    public :: engquist_osher_flux
    
contains
    
    ! Engquist-Osher flux implementation
    subroutine engquist_osher_flux(q_face_left, q_face_right, flux, config, mesh)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdConfigType), intent(in) :: config
        type(MeshType), intent(in) :: mesh
        
        integer :: i
        real(dp) :: c, cp, cm, u_L, u_R
        
        c = config%wave_speed
        
        do i = 1, mesh%nnodes
            cp = 0.5_dp * (c + abs(c))
            cm = 0.5_dp * (c - abs(c))
            u_L = q_face_left(i)
            u_R = q_face_right(i)
            flux(i) = cp * u_L + cm * u_R
        end do
    end subroutine engquist_osher_flux
    
end module engquist_osher_flux_module