! tests/test_registry.f90
program test_registry
    use registry_module
    implicit none
    
    type(component_info) :: info
    integer :: initial_size
    
    print *, "=== Fortran注册系统测试 ==="
    print *, ""
    
    ! 测试1: 基本注册
    print *, "测试1: 基本注册功能"
    print *, "---------------------"
    
    call component_registry%list_all()
    initial_size = component_registry%size()
    print *, "初始大小: ", initial_size
    
    ! 注册一些组件
    call register_component("reconstructor", "eno")
    call register_component("reconstructor", "weno3")
    call register_component("reconstructor", "weno5")
    call register_component("flux", "rusanov")
    call register_component("flux", "engquist-osher")
    call register_component("integrator", "rk1")
    call register_component("integrator", "rk2")
    call register_component("integrator", "rk3")
    
    print *, "注册后大小: ", component_registry%size()
    call component_registry%list_all()
    print *, ""
    
    ! 测试2: 重复注册
    print *, "测试2: 重复注册（应该显示警告）"
    print *, "-------------------------------"
    call register_component("reconstructor", "eno")
    print *, ""
    
    ! 测试3: 获取组件
    print *, "测试3: 获取组件功能"
    print *, "-------------------"
    
    ! 测试获取存在的组件
    info = component_registry%get("reconstructor", "weno3")
    print *, "获取 weno3: "
    call info%print()
    
    ! 测试获取不存在的组件
    info = component_registry%get("reconstructor", "non_existent")
    print *, "获取 non_existent (应该为空): "
    call info%print()
    print *, ""
    
    ! 测试4: 清空功能
    print *, "测试4: 清空注册表"
    print *, "-----------------"
    call component_registry%clear()
    print *, "清空后大小: ", component_registry%size()
    call component_registry%list_all()
    
    print *, ""
    print *, "=== 所有测试完成 ==="
    
end program test_registry