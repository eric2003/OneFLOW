! src/numerics/time_integration/rk3.f90
module rk3_integrator_module
    use base_modules, only: wp, ip
    use time_integration_base_module, only: time_integrator_base
    use solution_module, only: update_old_field
    
    implicit none
    private
    
    type, extends(time_integrator_base), public :: rk3_integrator
    contains
        procedure :: step => rk3_step
    end type rk3_integrator
    
contains

    subroutine rk3_step(this, dt)
        class(rk3_integrator), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        integer(ip) :: i, idx
        real(wp), allocatable :: u1(:), u2(:)
        integer(ip) :: n_cells
        
        ! 获取物理区域大小
        n_cells = this%domain%ied - this%domain%ist
        
        ! 分配中间解数组
        allocate(u1(this%domain%ist:this%domain%ied - 1))
        allocate(u2(this%domain%ist:this%domain%ied - 1))
        
        ! === 阶段1 ===
        call this%compute_residual()
        do i = this%domain%ist, this%domain%ied - 1
            idx = this%map_idx(i)
            u1(i) = this%solution%u(i) + dt * this%solution%res(idx)
        end do
        
        ! 应用边界条件到u1
        if (associated(this%bc)) then
            call this%bc%apply(u1, &
                               this%domain%nghosts, &
                               this%domain%ist, &
                               this%domain%ied - 1)
        end if
        
        ! === 阶段2 ===
        ! 将u1复制回当前解
        this%solution%u(this%domain%ist:this%domain%ied - 1) = &
            u1(this%domain%ist:this%domain%ied - 1)
        
        call this%compute_residual()
        
        do i = this%domain%ist, this%domain%ied - 1
            idx = this%map_idx(i)
            u2(i) = 0.75_wp * this%solution%un(i) + &
                    0.25_wp * this%solution%u(i) + &
                    0.25_wp * dt * this%solution%res(idx)
        end do
        
        ! 应用边界条件到u2
        if (associated(this%bc)) then
            call this%bc%apply(u2, &
                               this%domain%nghosts, &
                               this%domain%ist, &
                               this%domain%ied - 1)
        end if
        
        ! === 阶段3 ===
        ! 将u2复制回当前解
        this%solution%u(this%domain%ist:this%domain%ied - 1) = &
            u2(this%domain%ist:this%domain%ied - 1)
        
        call this%compute_residual()
        
        ! TVD RK3系数
        do i = this%domain%ist, this%domain%ied - 1
            idx = this%map_idx(i)
            this%solution%u(i) = (1.0_wp/3.0_wp) * this%solution%un(i) + &
                                 (2.0_wp/3.0_wp) * this%solution%u(i) + &
                                 (2.0_wp/3.0_wp) * dt * this%solution%res(idx)
        end do
        
        ! 最终边界条件
        call this%apply_boundary()
        
        ! 更新旧场
        call this%solution%update_old_field()
        
        ! 清理
        deallocate(u1, u2)
        
        if (this%config%verbose .and. mod(this%solution%current_step, 100) == 0) then
            print *, "  [RK3] Step completed"
        end if
    end subroutine rk3_step

end module rk3_integrator_module