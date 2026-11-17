! tests/test_factory_system.f90
program test_factory_system
    use registry_advanced_module
    use eno_reconstructor_module
    use rusanov_flux_module
    implicit none
    
    class(*), allocatable :: instance1, instance2
    type(eno_reconstructor), pointer :: eno_ptr
    type(rusanov_flux), pointer :: flux_ptr
    real(wp), allocatable :: q(:), qL(:), qR(:), flux(:)
    integer :: i
    
    print *, "=== Factory System Test ==="
    print *, ""
    
    ! 注册工厂
    print *, "1. Registering factories..."
    call register_factory_with_order("reconstructor", "eno", create_eno, 3)
    call register_factory("flux", "rusanov", create_rusanov)
    
    call component_registry%list_all()
    print *, ""
    
    ! 创建实例
    print *, "2. Creating instances..."
    call component_registry%create("reconstructor", "eno", instance1)
    call component_registry%create("flux", "rusanov", instance2)
    print *, ""
    
    ! 类型转换和测试
    print *, "3. Testing instances..."
    
    ! ENO重构器
    select type (inst => instance1)
    type is (eno_reconstructor)
        eno_ptr => inst
        call eno_ptr%info()
        
        ! 创建测试数据
        allocate(q(6), qL(5), qR(5))
        do i = 1, 6
            q(i) = real(i, wp)
        end do
        
        ! 计算
        call eno_ptr%compute(q, qL, qR)
        print *, "  q: ", q
        print *, "  qL: ", qL
        print *, "  qR: ", qR
        
        deallocate(q, qL, qR)
    class default
        print *, "[ERROR] Wrong type for eno_reconstructor"
    end select
    
    print *, ""
    
    ! Rusanov通量
    select type (inst => instance2)
    type is (rusanov_flux)
        flux_ptr => inst
        call flux_ptr%info()
        
        ! 创建测试数据
        allocate(qL(5), qR(5), flux(5))
        do i = 1, 5
            qL(i) = real(i, wp)
            qR(i) = real(i + 0.5, wp)
        end do
        
        ! 计算
        call flux_ptr%compute(qL, qR, flux)
        print *, "  qL: ", qL
        print *, "  qR: ", qR
        print *, "  flux: ", flux
        
        deallocate(qL, qR, flux)
    class default
        print *, "[ERROR] Wrong type for rusanov_flux"
    end select
    
    print *, ""
    print *, "=== Factory System Test Complete ==="
    
end program test_factory_system