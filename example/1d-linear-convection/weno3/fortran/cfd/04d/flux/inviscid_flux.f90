! ==================== inviscid_flux.f90 ====================
! Inviscid flux selector for OneFLOW-CFD
!
! Purpose: Select appropriate inviscid flux scheme based on configuration
! Author: OneFLOW-CFD Team
! Date: Created
! ==========================================================

module inviscid_flux_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    
    ! 导入具体通量模块
    use rusanov_flux_module, only: rusanov_flux
    use engquist_osher_flux_module, only: engquist_osher_flux
    
    implicit none
    
    private
    public :: inviscid_flux
    
contains
    
    ! Inviscid flux selection
    subroutine inviscid_flux(q_face_left, q_face_right, flux, config, mesh)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdConfigType), intent(in) :: config
        type(MeshType), intent(in) :: mesh
        
        ! 根据配置选择通量格式
        select case(config%flux_type)
        case(0)
            call rusanov_flux(q_face_left, q_face_right, flux, config, mesh)
        case(1)
            call engquist_osher_flux(q_face_left, q_face_right, flux, config, mesh)
        case default
            ! 默认使用Rusanov格式
            call rusanov_flux(q_face_left, q_face_right, flux, config, mesh)
        end select
    end subroutine inviscid_flux
    
end module inviscid_flux_module