! tests/test_physics_equations.f90
program test_physics_equations
  use precision, only: dp
  use linear_convection_equation, only: linear_convection_flux, linear_convection_speed
  use linear_convection_problem, only: linear_convection_prob, create_linear_convection_problem
  implicit none
  
  real(dp) :: u, f, a
  type(linear_convection_prob) :: prob
  real(dp), allocatable :: x(:), u_ic(:)
  integer :: i, nx = 10
  
  write(*,*) "=== 测试物理方程模块 ==="
  
  ! 测试线性对流方程
  write(*,*) "1. 测试线性对流方程..."
  u = 2.5_dp
  f = linear_convection_flux(u)
  a = linear_convection_speed()
  
  write(*,*) "   u = ", u
  write(*,*) "   F(u) = a*u = ", f
  write(*,*) "   a = ", a
  
  if (abs(f - 2.5_dp) < 1e-10_dp) then
    write(*,*) "   ✓ 通量计算正确"
  else
    write(*,*) "   ✗ 通量计算错误"
  end if
  
  ! 测试线性对流问题
  write(*,*) "2. 测试线性对流问题..."
  prob = create_linear_convection_problem(wave_speed=2.0_dp, &
                                          ic_type="step", &
                                          boundary_type="periodic")
  
  write(*,*) "   波速: ", prob%wave_speed
  write(*,*) "   初始条件类型: ", trim(prob%ic_type)
  write(*,*) "   边界类型: ", trim(prob%boundary_type)
  
  ! 测试通量方法
  u = 1.5_dp
  f = prob%flux(u)
  a = prob%speed()
  
  write(*,*) "   问题通量 F(u) = ", f, " (期望: 3.0)"
  write(*,*) "   问题波速 a = ", a, " (期望: 2.0)"
  
  if (abs(f - 3.0_dp) < 1e-10_dp) then
    write(*,*) "   ✓ 问题通量计算正确"
  else
    write(*,*) "   ✗ 问题通量计算错误"
  end if
  
  ! 测试初始条件
  write(*,*) "3. 测试初始条件..."
  allocate(x(nx), u_ic(nx))
  do i = 1, nx
    x(i) = 0.0_dp + (i-1) * 0.2_dp
  end do
  
  call prob%initial_condition(x, u_ic, nx)
  
  write(*,*) "   x 范围: ", x(1), " 到 ", x(nx)
  write(*,*) "   u_ic 范围: ", minval(u_ic), " 到 ", maxval(u_ic)
  
  ! 检查阶跃函数
  if (prob%ic_type == "step") then
    if (abs(u_ic(1) - 1.0_dp) < 1e-10_dp .and. &
        abs(u_ic(6) - 2.0_dp) < 1e-10_dp) then
      write(*,*) "   ✓ 阶跃初始条件正确"
    else
      write(*,*) "   ✗ 阶跃初始条件错误"
    end if
  end if
  
  deallocate(x, u_ic)
  
  write(*,*) "=== 物理方程模块测试完成 ==="
  write(*,*) "下一步: 实现初始条件工厂"
  
end program test_physics_equations