! ==================== flux.f90 ==========================
! Main flux computation module for OneFLOW-CFD
!
! Purpose: Provide unified interface for flux computations
! Author: OneFLOW-CFD Team
! Date: Created
! ========================================================

module flux_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    
    ! 重新导出所有公共过程
    use rusanov_flux_module, only: rusanov_flux
    use engquist_osher_flux_module, only: engquist_osher_flux
    use inviscid_flux_module, only: inviscid_flux
    
    implicit none
    
    ! 显式声明public接口，保持与原代码一致
    public :: rusanov_flux, engquist_osher_flux, inviscid_flux
    
end module flux_module