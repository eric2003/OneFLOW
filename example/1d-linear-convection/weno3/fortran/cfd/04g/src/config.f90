! ==================== config.f90 ====================
! Configuration module for OneFLOW-CFD

module config_module
    use kinds, only: dp
    implicit none
    
    private
    
    ! ===================================================================
    ! Configuration Type
    ! ===================================================================
    type, public :: CfdConfigType
        character(len=10) :: ic_type = "step"
        !character(len=10) :: ic_type = "sine"
		character(len=10) :: recon_scheme = "eno"
        integer :: flux_type = 0
        integer :: rk_order = 1
        integer :: spatial_order = 3
        real(dp) :: wave_speed = 1.0_dp
        real(dp) :: final_time = 0.625_dp
        real(dp) :: dt = 0.025_dp
        real(dp) :: cfl = 0.5_dp
    end type CfdConfigType
    
contains
    
end module config_module