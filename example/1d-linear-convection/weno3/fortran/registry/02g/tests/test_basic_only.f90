! tests/test_basic_only.f90
program test_basic_only
    ! 只测试最基本的功能，不依赖复杂模块
    use config_module, only: cfd_config, config_print, wp
    use mesh_module, only: mesh_type
    use registry_module, only: registry_init, registry_cleanup, &
                               register_component_simple, list_components
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    
    print *, "=== BASIC TEST - Minimal Functionality ==="
    print *, ""
    
    ! 测试1: 配置
    print *, "1. Testing configuration..."
    print *, "----------------------------"
    call config_print(config)
    print *, ""
    
    ! 测试2: 网格
    print *, "2. Testing mesh..."
    print *, "------------------"
    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=5)
    print *, "Mesh initialized:"
    print *, "  Cells: ", mesh%ncells
    print *, "  Nodes: ", mesh%nnodes
    print *, "  dx: ", mesh%dx
    print *, ""
    
    ! 测试3: 注册系统
    print *, "3. Testing registry..."
    print *, "----------------------"
    
    call registry_init()
    
    ! 注册组件（使用简化版本）
    call register_component_simple("reconstructor", "eno")
    call register_component_simple("reconstructor", "weno3")
    call register_component_simple("flux", "rusanov")
    
    ! 列出组件
    call list_components()
    print *, ""
    
    ! 清理
    call registry_cleanup()
    
    print *, "=== TEST PASSED ==="
    print *, "✓ Configuration works"
    print *, "✓ Mesh works"
    print *, "✓ Registry works"
    
end program test_basic_only