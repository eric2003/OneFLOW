! tests/test_factory_simple.f90 (修复版)
program test_factory_simple
    use base_modules, only: wp  ! ← 添加这行
    use config_module, only: cfd_config, config_print
    use mesh_module, only: mesh_type
    use reconstructor_base_module, only: reconstructor_base
    use eno_reconstructor_module, only: eno_reconstructor
    use weno3_reconstructor_module, only: weno3_reconstructor
    use flux_base_module, only: flux_calculator_base
    use rusanov_flux_module, only: rusanov_flux
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(eno_reconstructor) :: eno
    type(weno3_reconstructor) :: weno3
    type(rusanov_flux) :: rusanov
    
    print *, "=== Factory Pattern Simple Test ==="
    print *, ""
    
    ! Test 1: Basic systems
    print *, "1. Testing basic systems..."
    print *, "-----------------------------"
    call config_print(config)
    
    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=10)
    call mesh%print_info()
    print *, ""
    
    ! Test 2: Creating reconstructors
    print *, "2. Testing reconstructors..."
    print *, "------------------------------"
    
    ! 创建并测试ENO重构器
    print *, "Creating ENO reconstructor..."
    eno = eno_reconstructor()  ! 使用构造函数
    call eno%info()            ! 必须调用info方法
    
    print *, ""
    print *, "Creating WENO3 reconstructor..."
    weno3 = weno3_reconstructor()  ! 使用构造函数
    call weno3%info()              ! 必须调用info方法
    print *, ""
    
    ! Test 3: Creating flux calculator
    print *, "3. Testing flux calculator..."
    print *, "-------------------------------"
    
    print *, "Creating Rusanov flux calculator..."
    rusanov = rusanov_flux()  ! 使用构造函数
    call rusanov%info()       ! 必须调用info方法
    print *, ""
    
    print *, "=== Factory pattern simple test completed successfully ==="
    
end program test_factory_simple