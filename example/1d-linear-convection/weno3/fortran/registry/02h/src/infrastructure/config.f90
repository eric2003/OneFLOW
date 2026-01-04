! src/infrastructure/config.f90 (修复版)
module config_module
    use base_modules, only: wp, ip, max_name_len, cfd_config_base
    
    implicit none
    public :: wp, ip, cfd_config, config_print, config_with_reconstruction
    
    ! 扩展配置类型 - 添加物理相关字段
    type, extends(cfd_config_base) :: cfd_config
        ! 物理参数
        real(wp) :: left_boundary_value = 1.0_wp
        real(wp) :: right_boundary_value = 2.0_wp
        real(wp) :: domain_length = 2.0_wp
        
        ! 新增：物理模块相关配置
        real(wp) :: pulse_center = 0.5_wp      ! 高斯脉冲中心
        real(wp) :: pulse_width = 0.1_wp       ! 高斯脉冲宽度
        logical  :: enable_physics = .true.    ! 是否启用物理模块
    contains
        ! 新增：物理相关配置方法
        procedure :: set_physics_parameters
        procedure :: get_physics_info
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
        
        ! 新增：物理配置信息
        print *, "--- Physics Configuration ---"
        print *, "Equation type: ", trim(cfg%equation_type)
        print *, "Problem type: ", trim(cfg%problem_type)
        print *, "Domain length: ", cfg%domain_length
        print *, "Physics enabled: ", cfg%enable_physics
        
        if (cfg%ic_type == "gaussian") then
            print *, "Pulse center: ", cfg%pulse_center
            print *, "Pulse width: ", cfg%pulse_width
        end if
        
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
    
    ! ========== 新增：物理参数设置方法 ==========
    
    subroutine set_physics_parameters(this, equation_type, problem_type, &
                                      domain_length, enable_physics)
        class(cfd_config), intent(inout) :: this
        character(len=*), intent(in), optional :: equation_type, problem_type
        real(wp), intent(in), optional :: domain_length
        logical, intent(in), optional :: enable_physics
        
        if (present(equation_type)) then
            this%equation_type = trim(equation_type)
            if (this%verbose) then
                print *, "[CONFIG] Set equation type: ", trim(this%equation_type)
            end if
        end if
        
        if (present(problem_type)) then
            this%problem_type = trim(problem_type)
            if (this%verbose) then
                print *, "[CONFIG] Set problem type: ", trim(this%problem_type)
            end if
        end if
        
        if (present(domain_length)) then
            this%domain_length = domain_length
            if (this%verbose) then
                print *, "[CONFIG] Set domain length: ", this%domain_length
            end if
        end if
        
        if (present(enable_physics)) then
            this%enable_physics = enable_physics
            if (this%verbose) then
                print *, "[CONFIG] Physics module enabled: ", this%enable_physics
            end if
        end if
    end subroutine set_physics_parameters
    
    subroutine get_physics_info(this)
        class(cfd_config), intent(in) :: this
        
        print *, "=== Physics Configuration Info ==="
        print *, "Equation type: ", trim(this%equation_type)
        print *, "Problem type: ", trim(this%problem_type)
        print *, "Domain length: ", this%domain_length
        print *, "Wave speed: ", this%wave_speed
        print *, "Physics enabled: ", this%enable_physics
        
        if (this%ic_type == "gaussian") then
            print *, "Pulse parameters:"
            print *, "  Center: ", this%pulse_center
            print *, "  Width: ", this%pulse_width
        end if
        
        print *, "=================================="
    end subroutine get_physics_info

end module config_module