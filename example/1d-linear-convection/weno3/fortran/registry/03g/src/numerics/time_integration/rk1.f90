! src/numerics/time_integration/rk1.f90
module rk1_integrator_module
    use base_modules, only: wp, ip
    use time_integration_base_module, only: time_integrator_base
    use solution_module, only: update_old_field
    
    implicit none
    private
    
    type, extends(time_integrator_base), public :: rk1_integrator
    contains
        procedure :: step => rk1_step
    end type rk1_integrator
    
contains

    subroutine rk1_step(this, dt)
        class(rk1_integrator), intent(inout) :: this
        real(wp), intent(in) :: dt
        
        integer(ip) :: i, idx
        
        ! 计算残差
        call this%compute_residual()
        
        ! 更新解 (一阶RK = 欧拉方法)
        do i = this%domain%ist, this%domain%ied - 1
            idx = this%map_idx(i)
            this%solution%u(i) = this%solution%u(i) + dt * this%solution%res(idx)
        end do
        
        ! 应用边界条件
        call this%apply_boundary()
        
        ! 更新旧场
        call this%solution%update_old_field()
        
        if (this%config%verbose .and. mod(this%solution%current_step, 100) == 0) then
            print *, "  [RK1] Step completed"
        end if
    end subroutine rk1_step

end module rk1_integrator_module