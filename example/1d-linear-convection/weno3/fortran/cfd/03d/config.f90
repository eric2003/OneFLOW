! ==================== config.f90 ====================
! Configuration module for OneFLOW-CFD

module config_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    private
    
    ! ===================================================================
    ! Configuration Type
    ! ===================================================================
    type, public :: CfdConfigType
        character(len=10) :: recon_scheme = "eno"
        integer :: flux_type = 0
        integer :: rk_order = 1
        integer :: spatial_order = 3
        real(dp) :: wave_speed = 1.0_dp
        real(dp) :: final_time = 0.625_dp
        real(dp) :: dt = 0.025_dp
        real(dp) :: cfl = 0.5_dp
    end type CfdConfigType
    
    !public :: create_eno_config, create_weno_config
    !public :: create_weno_config
    
contains
    
    ! ===================================================================
    ! Configuration Creation Functions
    ! ===================================================================
    
    !!Create ENO configuration
    !function create_eno_config() result(config)
    !    type(CfdConfigType) :: config
    !    
    !    config%recon_scheme = "eno"
    !    config%spatial_order = 3
    !    config%flux_type = 0
    !    config%rk_order = 2
    !    config%wave_speed = 1.0_dp
    !    config%final_time = 0.625_dp
    !    config%cfl = 1.0_dp
    !    config%dt = 0.0025_dp
    !end function create_eno_config
    
    !! Create WENO configuration
    !function create_weno_config() result(config)
    !    type(CfdConfigType) :: config
    !    
    !    config%recon_scheme = "weno"
    !    config%spatial_order = 3
    !    config%flux_type = 0
    !    config%rk_order = 2
    !    config%wave_speed = 1.0_dp
    !    config%final_time = 0.625_dp
    !    config%cfl = 1.0_dp
    !    config%dt = 0.0025_dp
    !end function create_weno_config
    
end module config_module