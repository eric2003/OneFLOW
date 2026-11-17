! tests/test_component_manager.f90
program test_component_manager
    use, intrinsic :: iso_fortran_env, only: real64
    use config_module, only: cfd_config, config_print, config_with_reconstruction
    use component_manager_module, only: create_reconstructor, create_flux_calculator
    use component_manager_module, only: component_manager_info, validate_config
    use reconstructor_base_module, only: reconstructor_base
    use flux_base_module, only: flux_calculator_base
    implicit none
    
    type(cfd_config) :: config
    class(reconstructor_base), allocatable :: recon
    class(flux_calculator_base), allocatable :: flux
    integer :: status
    logical :: is_valid
    
    print *, "=== Component Manager Test ==="
    print *, ""
    
    ! 显示组件管理器信息
    call component_manager_info()
    print *, ""
    
    ! 测试1: 基本配置
    print *, "1. Testing basic ENO3 + Rusanov configuration..."
    print *, "-----------------------------------------------"
    
    config%verbose = .true.
    call config_print(config)
    
    ! 配置ENO3重构
    call config_with_reconstruction(config, "eno", 3)
    config%flux_type = "rusanov"
    config%wave_speed = 1.5_real64
    
    call config_print(config)
    print *, ""
    
    ! 验证配置
    is_valid = validate_config(config)
    if (is_valid) then
        print *, "[OK] Configuration is valid"
    else
        print *, "[ERROR] Configuration is invalid"
    end if
    print *, ""
    
    ! 测试2: 创建组件
    print *, "2. Testing component creation..."
    print *, "--------------------------------"
    
    ! 创建重构器（带状态检查）
    recon = create_reconstructor(config, status)
    if (status == 0) then
        print *, "[OK] Reconstructor created successfully"
        call recon%info()
    else
        print *, "[ERROR] Failed to create reconstructor, code:", status
    end if
    print *, ""
    
    ! 创建通量计算器
    flux = create_flux_calculator(config, status)
    if (status == 0) then
        print *, "[OK] Flux calculator created successfully"
        call flux%info()
    else
        print *, "[ERROR] Failed to create flux calculator, code:", status
    end if
    print *, ""
    
    ! 测试3: WENO3重构测试
    print *, "3. Testing WENO3 configuration..."
    print *, "---------------------------------"
    
    call config_with_reconstruction(config, "weno3", 3)
    
    is_valid = validate_config(config)
    if (is_valid) then
        print *, "[OK] WENO3 configuration is valid"
        
        ! 创建WENO3重构器
        recon = create_reconstructor(config)
        call recon%info()
    else
        print *, "[ERROR] WENO3 configuration is invalid"
    end if
    print *, ""
    
    ! 测试4: 错误配置测试
    print *, "4. Testing invalid configuration..."
    print *, "-----------------------------------"
    
    config%recon_scheme = "unknown_scheme"
    config%flux_type = "unknown_flux"
    
    is_valid = validate_config(config)
    if (.not. is_valid) then
        print *, "[OK] Invalid configuration correctly rejected"
    else
        print *, "[ERROR] Invalid configuration should have been rejected"
    end if
    
    ! 清理
    if (allocated(recon)) deallocate(recon)
    if (allocated(flux)) deallocate(flux)
    
    print *, ""
    print *, "=== Component manager test completed successfully ==="
    
end program test_component_manager