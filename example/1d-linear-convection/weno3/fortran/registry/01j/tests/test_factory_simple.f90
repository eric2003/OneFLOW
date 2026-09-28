program test_factory_simple
    use, intrinsic :: iso_fortran_env, only: real64
    use registry_module, only: initialize_registry, cleanup_registry, &
                               register_component_simple, has_component
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
    call config_print(config)
    call mesh%init(xmin=0.0_real64, xmax=1.0_real64, ncells=10)
    call mesh%print_info()
    print *, ""
    
    ! Test 2: Reconstructors - 必须调用info方法
    print *, "2. Testing reconstructors..."
    print *, "Creating ENO reconstructor..."
    
    ! 创建对象
    eno = eno_reconstructor()
    
    ! 必须调用info方法，否则链接器会认为不需要这些符号
    call eno%info()
    
    print *, "Creating WENO3 reconstructor..."
    weno3 = weno3_reconstructor()
    call weno3%info()
    print *, ""
    
    ! Test 3: Flux calculator
    print *, "3. Testing flux calculator..."
    rusanov = rusanov_flux()
    call rusanov%info()
    print *, ""
    
    ! Test 4: Registry
    print *, "4. Testing registry..."
    call initialize_registry(verbose=.true.)
    
    call register_component_simple("reconstructor", "eno")
    call register_component_simple("reconstructor", "weno3")
    call register_component_simple("flux", "rusanov")
    
    if (has_component("reconstructor", "eno")) then
        print *, "[OK] ENO reconstructor registered successfully"
    end if
    
    if (has_component("flux", "rusanov")) then
        print *, "[OK] Rusanov flux registered successfully"
    end if
    
    print *, ""
    call cleanup_registry()
    
    print *, "=== Factory pattern simple test completed successfully ==="
    
end program test_factory_simple