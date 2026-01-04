! tests/test_factory_simple.f90
program test_factory_simple
    use registry_module
    use config_module
    use mesh_module
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
    integer :: i
    
    print *, "=== Factory Pattern Simple Test ==="
    print *, ""
    
    ! Test 1: Configuration and mesh
    print *, "1. Testing basic systems..."
    print *, "-----------------------------"
    call config_print(config)
    call mesh%init(xmin=0.0_wp, xmax=1.0_wp, ncells=10)
    call mesh%print_info()
    print *, ""
    
    ! Test 2: Creating reconstructors directly
    print *, "2. Creating reconstructors..."
    print *, "------------------------------"
    eno = eno_reconstructor(name="ENO", order=3, epsilon=1.0e-6_wp)
    weno3 = weno3_reconstructor(name="WENO3", order=3, epsilon=1.0e-6_wp)
    
    call eno%info()
    call weno3%info()
    print *, ""
    
    ! Test 3: Creating flux calculator directly
    print *, "3. Creating flux calculator..."
    print *, "-------------------------------"
    rusanov = rusanov_flux(name="Rusanov", wave_speed_default=1.0_wp)
    call rusanov%info()
    print *, ""
    
    ! Test 4: Registry integration
    print *, "4. Testing registry..."
    print *, "----------------------"
    call initialize_registry(verbose=.true.)
    
    ! Register with simple method
    call register_component_simple("reconstructor", "eno")
    call register_component_simple("reconstructor", "weno3")
    call register_component_simple("flux", "rusanov")
    
    ! Check if registered
    if (has_component("reconstructor", "eno")) then
        print *, "✓ ENO reconstructor registered successfully"
    else
        print *, "✗ ENO reconstructor registration failed"
    end if
    
    if (has_component("flux", "rusanov")) then
        print *, "✓ Rusanov flux registered successfully"
    else
        print *, "✗ Rusanov flux registration failed"
    end if
    
    print *, ""
    call component_registry%list_all()
    print *, ""
    
    ! Test getting available components
    print *, "5. Testing component listing..."
    print *, "--------------------------------"
    call test_component_listing()
    print *, ""
    
    ! Cleanup
    call cleanup_registry()
    
    print *, "=== Factory pattern simple test completed successfully ==="
    
contains
    
    subroutine test_component_listing()
        character(len=:), allocatable :: recon_names(:)
        integer, allocatable :: recon_orders(:)
        integer :: i
        
        print *, "Available reconstructors:"
        call get_available_components("reconstructor", recon_names, recon_orders)
        
        if (allocated(recon_names)) then
            do i = 1, size(recon_names)
                print *, "  - ", trim(recon_names(i))
            end do
            print *, "Total reconstructors: ", size(recon_names)
        else
            print *, "  (no reconstructors found)"
        end if
    end subroutine test_component_listing
    
end program test_factory_simple