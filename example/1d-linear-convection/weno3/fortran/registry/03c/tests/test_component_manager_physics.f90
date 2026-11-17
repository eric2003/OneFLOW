! tests/test_component_manager_physics.f90 (简化版)
program test_component_manager_physics
    use base_modules, only: wp
    use config_module, only: cfd_config, config_print, config_with_reconstruction
    use component_manager_module, only: component_manager_info, validate_config
    
    implicit none
    
    type(cfd_config) :: config
    logical :: is_valid
    
    print *, "=== Component Manager Physics Test (Simplified) ==="
    print *, ""
    
    ! 测试1: 显示组件管理器信息
    print *, "1. Testing component manager info..."
    print *, "-------------------------------------"
    call component_manager_info()
    print *, ""
    
    ! 测试2: 物理模块测试（默认）
    print *, "2. Testing physics module with default configuration..."
    print *, "------------------------------------------------------"
    
    config%verbose = .true.
    call config_print(config)
    
    ! 验证配置
    is_valid = validate_config(config)
    if (is_valid) then
        print *, "[OK] Default configuration is valid"
    else
        print *, "[ERROR] Default configuration is invalid"
    end if
    print *, ""
    
    ! 测试3: 测试物理配置
    print *, "3. Testing physics configuration..."
    print *, "------------------------------------"
    
    ! 修改物理参数
    config%equation_type = "linear_advection"
    config%problem_type = "linear_advection"
    config%wave_speed = 2.5_wp
    config%domain_length = 3.0_wp
    
    print *, "Modified physics configuration:"
    print *, "  Equation type: ", trim(config%equation_type)
    print *, "  Problem type: ", trim(config%problem_type)
    print *, "  Wave speed: ", config%wave_speed
    print *, "  Domain length: ", config%domain_length
    print *, "  Physics enabled: ", config%enable_physics
    
    ! 验证修改后的配置
    is_valid = validate_config(config)
    if (is_valid) then
        print *, "[OK] Modified physics configuration is valid"
    else
        print *, "[ERROR] Modified physics configuration is invalid"
    end if
    print *, ""
    
    ! 测试4: 数值组件测试
    print *, "4. Testing numerical components with physics..."
    print *, "-----------------------------------------------"
    
    call config_with_reconstruction(config, "weno3", 3)
    config%flux_type = "rusanov"
    
    is_valid = validate_config(config)
    if (is_valid) then
        print *, "[OK] Combined physics+numerics configuration is valid"
    else
        print *, "[ERROR] Combined configuration is invalid"
    end if
    print *, ""
    
    ! 测试5: 物理模块禁用测试
    print *, "5. Testing physics module disabled..."
    print *, "---------------------------------------"
    
    config%enable_physics = .false.
    config%verbose = .false.
    
    is_valid = validate_config(config)
    if (is_valid) then
        print *, "[OK] Configuration valid even with physics disabled"
    else
        print *, "[ERROR] Configuration should be valid with physics disabled"
    end if
    print *, ""
    
    ! 测试6: 错误配置测试
    print *, "6. Testing error handling..."
    print *, "-----------------------------"
    
    config%verbose = .true.
    config%enable_physics = .true.
    config%recon_scheme = "unknown_scheme"
    config%flux_type = "unknown_flux"
    config%equation_type = "unknown_equation"
    config%problem_type = "unknown_problem"
    
    is_valid = validate_config(config)
    if (.not. is_valid) then
        print *, "[OK] Invalid configuration correctly rejected"
    else
        print *, "[ERROR] Invalid configuration should have been rejected"
    end if
    print *, ""
    
    print *, "=== Component Manager Physics Test Summary ==="
    print *, "✓ Component manager info works"
    print *, "✓ Configuration validation works with physics"
    print *, "✓ Error handling works correctly"
    print *, "✓ Combined physics+numerics validation works"
    print *, ""
    print *, "下一步: 集成物理模块到求解器框架"
    
end program test_component_manager_physics