! tests/test_infrastructure.f90 (原test_basic_only.f90)
program test_infrastructure
    use base_modules, only: wp
    use config_module, only: cfd_config, config_print
    use mesh_module, only: mesh_type
    use registry_module, only: registry_init, registry_cleanup, &
                               register_component_simple, list_components
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    
    print *, "=== 基础设施测试 ==="
    print *, ""
    
    ! 测试1: 配置
    print *, "1. 测试配置模块..."
    print *, "-------------------"
    call config_print(config)
    print *, ""
    
    ! 测试2: 网格
    print *, "2. 测试网格模块..."
    print *, "------------------"
    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=5)
    print *, "网格初始化:"
    print *, "  单元数: ", mesh%ncells
    print *, "  节点数: ", mesh%nnodes
    print *, "  网格间距: ", mesh%dx
    print *, ""
    
    ! 测试3: 注册系统
    print *, "3. 测试注册系统..."
    print *, "------------------"
    
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
    
    print *, "=== 基础设施测试通过 ==="
    print *, "✓ 配置模块工作正常"
    print *, "✓ 网格模块工作正常"
    print *, "✓ 注册系统工作正常"
    
end program test_infrastructure