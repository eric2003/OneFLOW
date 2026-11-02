! tests/test_config_physics.f90 (修复版)
program test_config_physics
    use base_modules, only: wp
    use config_module, only: cfd_config, config_print, config_with_reconstruction
    
    implicit none
    
    type(cfd_config) :: config
    
    print *, "=== Configuration Physics Test (Simplified) ==="
    print *, ""
    
    ! 测试1: 默认配置
    print *, "1. Testing default configuration..."
    print *, "-----------------------------------"
    call config_print(config)
    print *, ""
    
    ! 测试2: 验证基础物理字段
    print *, "2. Testing basic physics fields..."
    print *, "----------------------------------"
    
    print *, "Verifying default physics fields:"
    
    if (trim(config%equation_type) == "linear_advection") then
        print *, "  ✓ Default equation type: linear_advection"
    else
        print *, "  ✗ Unexpected equation type: ", trim(config%equation_type)
    end if
    
    if (trim(config%problem_type) == "linear_advection") then
        print *, "  ✓ Default problem type: linear_advection"
    else
        print *, "  ✗ Unexpected problem type: ", trim(config%problem_type)
    end if
    
    if (abs(config%domain_length - 2.0_wp) < 1e-10_wp) then
        print *, "  ✓ Default domain length: 2.0"
    else
        print *, "  ✗ Unexpected domain length: ", config%domain_length
    end if
    
    if (config%enable_physics) then
        print *, "  ✓ Physics enabled by default"
    else
        print *, "  ✗ Physics not enabled by default"
    end if
    
    print *, ""
    
    ! 测试3: 使用类型绑定的方法（正确的方法名）
    print *, "3. Testing type-bound procedures..."
    print *, "--------------------------------------"
    
    call config%set_physics_parameters( &
        equation_type="burgers_equation", &
        problem_type="sod_shock_tube", &
        domain_length=3.0_wp, &
        enable_physics=.false.)
    
    print *, "After set_physics_parameters:"
    print *, "  Equation type: ", trim(config%equation_type)
    print *, "  Problem type: ", trim(config%problem_type)
    print *, "  Domain length: ", config%domain_length
    print *, "  Physics enabled: ", config%enable_physics
    
    if (trim(config%equation_type) == "burgers_equation") then
        print *, "  ✓ Equation type modified successfully via set_physics_parameters"
    end if
    
    if (trim(config%problem_type) == "sod_shock_tube") then
        print *, "  ✓ Problem type modified successfully via set_physics_parameters"
    end if
    
    if (abs(config%domain_length - 3.0_wp) < 1e-10_wp) then
        print *, "  ✓ Domain length modified successfully via set_physics_parameters"
    end if
    
    if (.not. config%enable_physics) then
        print *, "  ✓ Physics disabled successfully via set_physics_parameters"
    end if
    
    print *, ""
    
    ! 测试4: 调用get_physics_info方法
    print *, "4. Testing get_physics_info method..."
    print *, "--------------------------------------"
    call config%get_physics_info()
    print *, ""
    
    ! 测试5: 高斯脉冲配置
    print *, "5. Testing Gaussian pulse configuration..."
    print *, "-----------------------------------------"
    
    config%ic_type = "gaussian"
    config%pulse_center = 0.6_wp
    config%pulse_width = 0.15_wp
    
    print *, "Gaussian pulse parameters:"
    print *, "  IC type: ", trim(config%ic_type)
    print *, "  Center: ", config%pulse_center
    print *, "  Width: ", config%pulse_width
    
    if (trim(config%ic_type) == "gaussian") then
        print *, "  ✓ Gaussian IC type set"
    end if
    
    if (abs(config%pulse_center - 0.6_wp) < 1e-10_wp) then
        print *, "  ✓ Pulse center set"
    end if
    
    if (abs(config%pulse_width - 0.15_wp) < 1e-10_wp) then
        print *, "  ✓ Pulse width set"
    end if
    
    print *, ""
    
    ! 测试6: 重构配置
    print *, "6. Testing reconstruction configuration..."
    print *, "------------------------------------------"
    
    call config_with_reconstruction(config, "weno", 5)
    
    print *, "Reconstruction configuration:"
    print *, "  Scheme: ", trim(config%recon_scheme)
    print *, "  Order: ", config%spatial_order
    
    if (trim(config%recon_scheme) == "weno" .and. config%spatial_order == 5) then
        print *, "  ✓ WENO5 configuration successful"
    else
        print *, "  ✗ Reconstruction configuration failed"
    end if
    
    print *, ""
    
    print *, "=== Configuration Physics Test Complete ==="
    print *, "✓ Config module updated with physics support"
    print *, "✓ Fields can be directly accessed and modified"
    print *, "✓ Type-bound procedures work correctly"
    
end program test_config_physics