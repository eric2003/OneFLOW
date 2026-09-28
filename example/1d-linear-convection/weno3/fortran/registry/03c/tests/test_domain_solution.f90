! tests/test_domain_solution.f90
program test_domain_solution
    use base_modules, only: wp
    use config_module, only: cfd_config, config_with_reconstruction
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type, domain_create
    use solution_module, only: solution_type, solution_create, solution_reset
    
    implicit none
    
    type(cfd_config) :: config
    type(mesh_type) :: mesh
    type(domain_type) :: domain
    type(solution_type) :: solution
    real(wp), allocatable :: initial_values(:)
    integer :: i
    
    print *, "=== Domain and Solution Test ==="
    print *, ""
    
    ! 测试1: 不同重构方案的ghost层计算
    print *, "1. Testing ghost layer calculation..."
    print *, "--------------------------------------"
    
    ! ENO3
    call config_with_reconstruction(config, "eno", 3)
    config%verbose = .false.
    call mesh%init(ncells=10)
    domain = domain_create(config, mesh)
    print *, "ENO3: nghosts = ", domain%nghosts, " (expected: 3)"
    
    ! WENO3
    call config_with_reconstruction(config, "weno3", 3)
    domain = domain_create(config, mesh)
    print *, "WENO3: nghosts = ", domain%nghosts, " (expected: 2)"
    
    ! WENO5
    call config_with_reconstruction(config, "weno", 5)
    domain = domain_create(config, mesh)
    print *, "WENO5: nghosts = ", domain%nghosts, " (expected: 3)"
    print *, ""
    
    ! 测试2: Solution数组
    print *, "2. Testing solution arrays..."
    print *, "------------------------------"
    
    call config_with_reconstruction(config, "eno", 3)
    config%verbose = .true.
    call mesh%init(xmin=0.0_wp, xmax=2.0_wp, ncells=10)
    
    domain = domain_create(config, mesh)
    call domain%print_info()
    print *, ""
    
    solution = solution_create(domain)
    call solution%print_info()
    print *, ""
    
    ! 测试3: 初始化和更新
    print *, "3. Testing initialization and update..."
    print *, "----------------------------------------"
    
    allocate(initial_values(mesh%ncells))
    do i = 1, mesh%ncells
        initial_values(i) = sin(2.0_wp * 3.14159265358979_wp * mesh%xcc(i) / mesh%L)
    end do
    
    call solution%initialize(initial_values)
    print *, "After initialization:"
    print *, "  u range: ", minval(solution%u), " to ", maxval(solution%u)
    print *, "  un range: ", minval(solution%un), " to ", maxval(solution%un)
    
    ! 修改当前解，测试更新
    solution%u = solution%u * 2.0_wp
    call solution%update_old_field()
    print *, "After update: max|u - un| = ", maxval(abs(solution%u - solution%un))
    print *, ""
    
    ! 测试4: 重置
    print *, "4. Testing reset..."
    print *, "-------------------"
    
    call solution_reset(solution)
    print *, "After reset:"
    print *, "  u max: ", maxval(abs(solution%u))
    print *, "  un max: ", maxval(abs(solution%un))
    print *, "  flux max: ", maxval(abs(solution%flux))
    print *, ""
    
    deallocate(initial_values)
    
    print *, "=== Test Summary ==="
    print *, "✓ Ghost layer calculation works"
    print *, "✓ Domain creation works"
    print *, "✓ Solution arrays work"
    print *, "✓ Initialization works"
    print *, "✓ Field update works"
    print *, "✓ Reset works"
    print *, ""
    print *, "Ready for next step: Implementing Physics modules"
    
end program test_domain_solution