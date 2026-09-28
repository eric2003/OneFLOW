! tests/test_factory.f90
program test_factory
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
    class(*), allocatable :: instance
    class(reconstructor_base), allocatable :: recon
    class(flux_calculator_base), allocatable :: flux_calc
    real(wp), allocatable :: q(:), qL(:), qR(:), flux(:)
    integer :: i, n

    print *, "=== Factory Pattern Test ==="
    print *, ""

    ! Initialize systems
    call initialize_registry(verbose=.true.)

    ! Register components with factories
    print *, "1. Registering components with factories..."
    call register_component_with_factory("reconstructor", "eno", create_eno, 3)
    call register_component_with_factory("reconstructor", "weno3", create_weno3, 3)
    call register_component_with_factory("flux", "rusanov", create_rusanov)

    call component_registry%list_all()
    print *, ""

    ! Test creating ENO reconstructor
    print *, "2. Creating ENO reconstructor..."
    call create_component("reconstructor", "eno", instance)

    if (allocated(instance)) then
        select type (inst => instance)
        type is (eno_reconstructor)
            allocate(recon, source=inst)
            print *, "ENO reconstructor created successfully"
            call recon%info()
            print *, ""

            ! Test reconstruction
            n = 10
            allocate(q(0:n+1), qL(n), qR(n))

            ! Initialize test data (sine wave)
            do i = 0, n+1
                q(i) = sin(2.0_wp * 3.141592653589793_wp * real(i-1, wp) / real(n, wp))
            end do

            print *, "Testing ENO reconstruction..."
            call recon%reconstruct(q, qL, qR)

            print *, "q (internal):"
            do i = 1, n
                write(*, '(I3, F10.6)') i, q(i)
            end do

            print *, "qL (left interface values):"
            do i = 1, n
                write(*, '(I3, F10.6)') i, qL(i)
            end do

            deallocate(q, qL, qR)
        class default
            print *, "[ERROR] Wrong type for ENO reconstructor"
        end select
    else
        print *, "[ERROR] Failed to create ENO reconstructor"
    end if
    print *, ""

    ! Test creating WENO3 reconstructor
    print *, "3. Creating WENO3 reconstructor..."
    if (allocated(instance)) deallocate(instance)
    call create_component("reconstructor", "weno3", instance)

    if (allocated(instance)) then
        select type (inst => instance)
        type is (weno3_reconstructor)
            if (allocated(recon)) deallocate(recon)
            allocate(recon, source=inst)
            print *, "WENO3 reconstructor created successfully"
            call recon%info()
        class default
            print *, "[ERROR] Wrong type for WENO3 reconstructor"
        end select
    else
        print *, "[ERROR] Failed to create WENO3 reconstructor"
    end if
    print *, ""

    ! Test creating Rusanov flux calculator
    print *, "4. Creating Rusanov flux calculator..."
    if (allocated(instance)) deallocate(instance)
    call create_component("flux", "rusanov", instance)

    if (allocated(instance)) then
        select type (inst => instance)
        type is (rusanov_flux)
            allocate(flux_calc, source=inst)
            print *, "Rusanov flux calculator created successfully"
            call flux_calc%info()
            print *, ""

            ! Test flux computation
            n = 5
            allocate(qL(n), qR(n), flux(n))

            ! Initialize test data
            do i = 1, n
                qL(i) = 1.0_wp + 0.1_wp * real(i-1, wp)
                qR(i) = 1.0_wp + 0.1_wp * real(i, wp)
            end do

            print *, "Testing Rusanov flux computation..."
            call flux_calc%compute(qL, qR, flux, 1.0_wp)

            print *, "qL:"
            do i = 1, n
                write(*, '(I3, F10.6)') i, qL(i)
            end do

            print *, "qR:"
            do i = 1, n
                write(*, '(I3, F10.6)') i, qR(i)
            end do

            print *, "Flux:"
            do i = 1, n
                write(*, '(I3, F10.6)') i, flux(i)
            end do

            deallocate(qL, qR, flux)
        class default
            print *, "[ERROR] Wrong type for Rusanov flux"
        end select
    else
        print *, "[ERROR] Failed to create Rusanov flux calculator"
    end if
    print *, ""

    ! Cleanup
    call cleanup_registry()

    print *, "=== Factory pattern test completed successfully ==="

end program test_factory