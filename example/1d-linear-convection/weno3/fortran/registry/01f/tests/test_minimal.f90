! tests/test_minimal.f90
program test_minimal
    use registry_module
    use config_module
    use mesh_module
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    integer :: i  ! Declare loop variable here
    
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
    
    call initialize_registry(verbose=.true.)
    
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
    
    call component_registry%list_all()
    print *, "Registry size: ", component_registry%size()
    print *, ""
    
    ! Test component lookup
    print *, "4. Testing component lookup"
    print *, "---------------------------"
    if (has_component("reconstructor", "eno")) then
        print *, "Found: reconstructor.eno"
    else
        print *, "Not found: reconstructor.eno"
    end if
    
    if (has_component("reconstructor", "unknown")) then
        print *, "Found: reconstructor.unknown"
    else
        print *, "Not found: reconstructor.unknown"
    end if
    print *, ""
    
    ! Test getting available components
    print *, "5. Testing available components"
    print *, "--------------------------------"
    
    ! Use a separate block to avoid allocation issues
    call test_available_components()
    print *, ""
    
    ! Cleanup
    call cleanup_registry()
    
    print *, "=== Minimal test completed successfully ==="
    
contains
    
    ! Internal subroutine for testing available components
    subroutine test_available_components()
        character(len=:), allocatable :: names(:)
        integer, allocatable :: orders(:)
        integer :: j
        
        call get_available_components("reconstructor", names, orders)
        
        if (allocated(names)) then
            print *, "Reconstructors:"
            do j = 1, size(names)
                print *, "  - ", trim(names(j))
                if (allocated(orders)) then
                    print *, "    Order: ", orders(j)
                end if
            end do
        else
            print *, "No reconstructors found"
        end if
    end subroutine test_available_components
    
end program test_minimal