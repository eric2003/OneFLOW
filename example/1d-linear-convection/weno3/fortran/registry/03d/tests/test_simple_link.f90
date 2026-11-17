! tests/test_simple_link.f90
program test_simple_link
    use base_modules, only: wp  ! ← 添加这行
    use registry_module
    use config_module
    use mesh_module
    implicit none

    type(cfd_config) :: config
    type(mesh_type) :: mesh
    integer :: i

    print *, "=== Minimal Functionality Test ==="
    print *, ""

    ! Test 1: Configuration system
    print *, "1. Testing configuration system"
    print *, "--------------------------------"

    call config_print(config)

    call config_with_reconstruction(config, "eno", 3)

    call config_print(config)
    print *, ""

    ! Test 2: Mesh system
    print *, "2. Testing mesh system"
    print *, "----------------------"

    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=10)
    call mesh%print_info()
    print *, ""

    ! Test 3: Registry system
    print *, "3. Testing registry system"
    print *, "--------------------------"

    call registry_init()

    ! Register some components
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
    print *, "Registry size: ", registry_get_size()
    print *, ""

    ! Test component lookup
    print *, "4. Testing component lookup"
    print *, "---------------------------"
    if (has_component_simple("reconstructor", "eno")) then
        print *, "Found: reconstructor.eno"
    else
        print *, "Not found: reconstructor.eno"
    end if

    if (has_component_simple("reconstructor", "unknown")) then
        print *, "Found: reconstructor.unknown"
    else
        print *, "Not found: reconstructor.unknown"
    end if
    print *, ""

    ! Cleanup
    call registry_cleanup()

    print *, "=== Minimal test completed successfully ==="

end program test_simple_link