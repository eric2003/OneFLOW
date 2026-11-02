! src/numerics/time_integration/rk2.f90
module rk2_integrator_module
    use base_modules, only: wp, ip
    use time_integration_base_module, only: time_integrator_base
    use solution_module, only: update_old_field
    
    implicit none
    private
    
    type, extends(time_integrator_base), public :: rk2_integrator
    contains
        procedure :: step => rk2_step
    end type rk2_integrator
    
contains

    subroutine rk2_step(this, dt)
        class(rk2_integrator), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        integer(ip) :: i, idx
        real(wp), allocatable :: u_pred(:)
        integer(ip) :: n_cells
        
        ! 获取物理区域大小
        n_cells = this%domain%ied - this%domain%ist
        
        ! 分配预测解数组
        allocate(u_pred(this%domain%ist:this%domain%ied - 1))
        
        ! 阶段1: 计算预测值
        call this%compute_residual()
        do i = this%domain%ist, this%domain%ied - 1
            idx = this%map_idx(i)
            u_pred(i) = this%solution%u(i) + dt * this%solution%res(idx)
        end do
        
        ! 应用边界条件到预测值
        if (associated(this%bc)) then
            call this%bc%apply(u_pred, &
                               this%domain%nghosts, &
                               this%domain%ist, &
                               this%domain%ied - 1)
        end if
        
        ! 阶段2: 计算修正值
        ! 将预测值复制回当前解
        this%solution%u(this%domain%ist:this%domain%ied - 1) = &
            u_pred(this%domain%ist:this%domain%ied - 1)
        
        ! 再次计算残差
        call this%compute_residual()
        
        ! 计算最终值 (Heun方法: u^{n+1} = 0.5*u^n + 0.5*(u_pred + dt*res(u_pred)))
        do i = this%domain%ist, this%domain%ied - 1
            idx = this%map_idx(i)
            this%solution%u(i) = 0.5_wp * this%solution%un(i) + &
                                 0.5_wp * this%solution%u(i) + &
                                 0.5_wp * dt * this%solution%res(idx)
        end do
        
        ! 应用边界条件
        call this%apply_boundary()
        
        ! 更新旧场
        call this%solution%update_old_field()
        
        ! 清理
        deallocate(u_pred)
        
        if (this%config%verbose .and. mod(this%solution%current_step, 100) == 0) then
            print *, "  [RK2] Step completed"
        end if
    end subroutine rk2_step

end module rk2_integrator_module