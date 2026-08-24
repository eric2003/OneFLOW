! ==================== flux.f90 ====================
! Flux computation module for OneFLOW-CFD

module flux_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    implicit none
    
    private
    
    ! ===================================================================
    ! Flux Computation Interface
    ! ===================================================================
    ! 这里不需要定义新类型，只需要提供通量计算函数
    
    public :: rusanov_flux, engquist_osher_flux, inviscid_flux
    
contains
    
    ! ===================================================================
    ! Flux Functions
    ! ===================================================================
    
    ! Rusanov flux
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
    
    ! Engquist-Osher flux
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
    
    ! Inviscid flux selection
    subroutine inviscid_flux(q_face_left, q_face_right, flux, config, mesh)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdConfigType), intent(in) :: config
        type(MeshType), intent(in) :: mesh
        
        if (config%flux_type == 0) then
            call rusanov_flux(q_face_left, q_face_right, flux, config, mesh)
        else
            call engquist_osher_flux(q_face_left, q_face_right, flux, config, mesh)
        end if
    end subroutine inviscid_flux
    
end module flux_module