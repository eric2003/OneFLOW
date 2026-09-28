! src/base/modules.f90
module base_modules
    use, intrinsic :: iso_fortran_env, only: real64, int32
    
    implicit none
    
    ! 公开所有基本类型
    public :: wp, ip, string_len, max_name_len
    public :: cfd_config_base, component_info
    
    ! 工作精度类型
    integer, parameter :: wp = real64
    integer, parameter :: ip = int32
    
    ! 字符串长度常量（使用数值常量）
    integer, parameter :: string_len = 100
    integer, parameter :: max_name_len = 32
    
    ! 基础配置类型
    type :: cfd_config_base
        character(len=max_name_len) :: ic_type = "step"
        character(len=max_name_len) :: recon_scheme = "eno"
        character(len=max_name_len) :: flux_type = "rusanov"
        integer(ip) :: rk_order = 1
        real(wp) :: wave_speed = 1.0_wp
        real(wp) :: final_time = 0.625_wp
        real(wp) :: dt = 0.025_wp
        character(len=max_name_len) :: boundary_type = "periodic"
        integer(ip) :: spatial_order = 2
        character(len=max_name_len) :: equation_type = "linear_advection"
        character(len=max_name_len) :: problem_type = "linear_advection"
        logical :: verbose = .true.
    end type cfd_config_base
    
    ! 组件信息类型
    type :: component_info
        character(len=max_name_len) :: category = ""
        character(len=max_name_len) :: name = ""
        integer(ip) :: order = 0
    end type component_info
    
end module base_modules