! tests/test_registry_basic.f90
program test_registry_basic
    use registry_module
    implicit none
    
    type(component_info) :: info
    logical :: found
    
    print *, "=== 注册系统基础测试 ==="
    print *, ""
    
    ! 1. 初始化
    print *, "1. 初始化注册表"
    print *, "----------------"
    call initialize_registry(verbose=.true.)
    print *, "初始大小: ", component_registry%size()
    print *, ""
    
    ! 2. 注册组件
    print *, "2. 注册组件"
    print *, "------------"
    call register_component_simple("reconstructor", "eno")
    call register_component_simple("reconstructor", "weno3")
    call register_component_simple("reconstructor", "weno5")
    call register_component_simple("flux", "rusanov")
    call register_component_simple("flux", "engquist-osher")
    call register_component_simple("boundary", "periodic")
    call register_component_simple("boundary", "dirichlet")
    call register_component_simple("integrator", "rk1")
    call register_component_simple("integrator", "rk2")
    call register_component_simple("integrator", "rk3")
    
    call component_registry%list_all()
    print *, "注册后大小: ", component_registry%size()
    print *, ""
    
    ! 3. 查询组件
    print *, "3. 查询组件"
    print *, "------------"
    
    ! 测试存在的组件
    found = has_component("reconstructor", "eno")
    print *, "has_component('reconstructor', 'eno') = ", found
    
    info = component_registry%get("reconstructor", "eno")
    print *, "get('reconstructor', 'eno'):"
    call info%print()
    
    ! 测试不存在的组件
    found = has_component("reconstructor", "unknown")
    print *, "has_component('reconstructor', 'unknown') = ", found
    
    info = component_registry%get("reconstructor", "unknown")
    print *, "get('reconstructor', 'unknown'):"
    call info%print()
    print *, ""
    
    ! 4. 清理
    print *, "4. 清理注册表"
    print *, "--------------"
    call cleanup_registry()
    print *, "清理后大小: ", component_registry%size()
    call component_registry%list_all()
    
    print *, ""
    print *, "=== 基础测试完成 ==="
    
end program test_registry_basic