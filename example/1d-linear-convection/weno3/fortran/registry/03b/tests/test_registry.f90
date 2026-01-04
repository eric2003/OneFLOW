! tests/test_registry.f90 (原test_minimal_simple.f90)
program test_registry
    use base_modules, only: wp
    use registry_module
    use config_module
    use mesh_module
    
    implicit none

    type(cfd_config) :: config
    type(mesh_type) :: mesh
    integer :: i

    print *, "=== 注册系统功能测试 ==="
    print *, ""

    ! 测试1: 配置系统
    print *, "1. 测试配置系统"
    print *, "--------------"

    call config_print(config)

    call config_with_reconstruction(config, "eno", 3)

    call config_print(config)
    print *, ""

    ! 测试2: 网格系统
    print *, "2. 测试网格系统"
    print *, "--------------"

    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=10)
    call mesh%print_info()
    print *, ""

    ! 测试3: 注册系统
    print *, "3. 测试注册系统"
    print *, "--------------"

    call registry_init()

    ! 注册组件
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

    call list_components()
    print *, "注册表大小: ", registry_get_size()
    print *, ""

    ! 测试组件查找
    print *, "4. 测试组件查找"
    print *, "--------------"
    
    if (has_component_simple("reconstructor", "eno")) then
        print *, "找到: reconstructor.eno"
    else
        print *, "未找到: reconstructor.eno"
    end if

    if (has_component_simple("reconstructor", "unknown")) then
        print *, "找到: reconstructor.unknown"
    else
        print *, "未找到: reconstructor.unknown"
    end if
    print *, ""

    ! 测试获取可用组件
    print *, "5. 测试注册系统功能"
    print *, "------------------"
    print *, "注册表已初始化: ", registry_is_initialized()
    print *, "组件数量: ", registry_get_size()
    print *, ""

    ! 清理
    call registry_cleanup()

    print *, "=== 注册系统测试完成 ==="

end program test_registry