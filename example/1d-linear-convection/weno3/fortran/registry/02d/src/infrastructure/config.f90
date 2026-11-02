! src/infrastructure/config.f90
module config_module
    use base_modules, only: wp, ip, max_name_len, cfd_config_base
    
    implicit none
    public :: wp, ip, cfd_config, config_print, config_with_reconstruction
    
    ! 扩展配置类型
    type, extends(cfd_config_base) :: cfd_config
        real(wp) :: left_boundary_value = 1.0_wp
        real(wp) :: right_boundary_value = 2.0_wp
        real(wp) :: domain_length = 2.0_wp
    end type cfd_config
    
contains

    subroutine config_print(cfg)
        type(cfd_config), intent(in) :: cfg
        
        print *, "=== CFD Configuration ==="
        print *, "Initial condition: ", trim(cfg%ic_type)
        print *, "Reconstruction: ", trim(cfg%recon_scheme), " (order:", cfg%spatial_order, ")"
        print *, "Flux type: ", trim(cfg%flux_type)
        print *, "Time integration: RK", cfg%rk_order
        print *, "Wave speed: ", cfg%wave_speed
        print *, "Final time: ", cfg%final_time
        print *, "Time step: ", cfg%dt
        print *, "Boundary: ", trim(cfg%boundary_type)
        print *, "==============================="
    end subroutine config_print

    subroutine config_with_reconstruction(cfg, scheme, order)
        type(cfd_config), intent(inout) :: cfg
        character(len=*), intent(in) :: scheme
        integer, optional, intent(in) :: order
        
        integer :: i
        
        ! 转换为小写
        cfg%recon_scheme = scheme
        do i = 1, len_trim(cfg%recon_scheme)
            if (cfg%recon_scheme(i:i) >= 'A' .and. cfg%recon_scheme(i:i) <= 'Z') then
                cfg%recon_scheme(i:i) = char(ichar(cfg%recon_scheme(i:i)) + 32)
            end if
        end do
        
        ! 设置阶数
        if (present(order)) then
            cfg%spatial_order = order
        else
            if (index(cfg%recon_scheme, 'weno') > 0) then
                cfg%spatial_order = 5
            else if (trim(cfg%recon_scheme) == 'eno') then
                cfg%spatial_order = 3
            end if
        end if
        
        if (cfg%verbose) then
            print *, "[CONFIG] Reconstruction: ", trim(cfg%recon_scheme), &
                     " Order: ", cfg%spatial_order
        end if
    end subroutine config_with_reconstruction

end module config_module