! tests/test_registry.f90
program test_registry
    use registry_module
    implicit none
    
    type(component_info) :: info
    
    print *, "=== Fortran Registry System Test ==="
    print *, ""
    
    ! 初始化注册表
    call initialize_registry(verbose=.true.)
    print *, ""
    
    ! 测试1: 基本注册
    print *, "Test 1: Basic Registration"
    print *, "--------------------------"
    
    call component_registry%list_all()
    print *, "Initial size: ", component_registry%size()
    
    ! 使用两种方式注册组件
    call register_component("reconstructor", "eno")
    call register_component("reconstructor", "weno3")  ! 简单版本
    call register_component("reconstructor", "weno5")
    call register_component("flux", "rusanov")
    call register_component("flux", "engquist-osher")
    
    print *, "After registration size: ", component_registry%size()
    call component_registry%list_all()
    print *, ""
    
    ! 测试2: 重复注册
    print *, "Test 2: Duplicate Registration"
    print *, "-------------------------------"
    call register_component("reconstructor", "eno")
    print *, ""
    
    ! 测试3: 获取组件
    print *, "Test 3: Component Retrieval"
    print *, "---------------------------"
    
    info = component_registry%get("reconstructor", "eno")
    print *, "Get eno: "
    call info%print()
    
    info = component_registry%get("reconstructor", "non_existent")
    print *, "Get non_existent (should be empty): "
    call info%print()
    print *, ""
    
    ! 测试4: 清空注册表
    print *, "Test 4: Clear Registry"
    print *, "----------------------"
    call component_registry%clear()
    print *, "Size after clear: ", component_registry%size()
    call component_registry%list_all()
    print *, ""
    
    ! 测试5: 重新初始化并注册
    print *, "Test 5: Re-initialize and register"
    print *, "----------------------------------"
    call initialize_registry(verbose=.false.)
    call register_component("integrator", "rk1")
    call register_component("integrator", "rk2")
    call register_component("integrator", "rk3")
    call component_registry%list_all()
    print *, ""
    
    ! 清理
    call cleanup_registry()
    
    print *, "=== All Tests Completed Successfully ==="
    
end program test_registry