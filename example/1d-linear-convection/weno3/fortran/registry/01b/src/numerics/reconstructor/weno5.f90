! src/numerics/reconstructor/weno5.f90
module weno5_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use reconstructor_base_module, only: reconstructor_base, wp
    use registry_module, only: register_component_with_factory
    implicit none
    
    private
    public :: wp, weno5_reconstructor, create_weno5
    
    type, extends(reconstructor_base) :: weno5_reconstructor
    contains
        procedure :: reconstruct => weno5_reconstruct
        procedure :: info => weno5_info
    end type weno5_reconstructor
    
contains
    
    ! 工厂函数
    subroutine create_weno5(instance)
        class(*), allocatable, intent(out) :: instance
        type(weno5_reconstructor), allocatable :: weno5
        
        allocate(weno5)
        weno5%name = "WENO5"
        weno5%order = 5
        weno5%epsilon = 1.0e-6_wp
        
        call move_alloc(weno5, instance)
    end subroutine create_weno5
    
    ! WENO-5重构实现（简化版本）
    subroutine weno5_reconstruct(this, q, qL, qR)
        class(weno5_reconstructor), intent(in) :: this
        real(wp), intent(in) :: q(:)
        real(wp), intent(out) :: qL(:), qR(:)
        
        integer :: i, n
        real(wp) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, alpha
        real(wp) :: w0, w1, w2, q0, q1, q2
        real(wp) :: v0, v1, v2, v3, v4
        
        n = size(q) - 4  ! 需要4个ghost cells
        
        do i = 3, n+2
            ! 获取模板值
            v0 = q(i-2)
            v1 = q(i-1)
            v2 = q(i)
            v3 = q(i+1)
            v4 = q(i+2)
            
            ! 计算左界面值 qL(i-2)
            beta0 = (13.0_wp/12.0_wp)*(v0 - 2.0_wp*v1 + v2)**2 &
                  + (1.0_wp/4.0_wp)*(v0 - 4.0_wp*v1 + 3.0_wp*v2)**2
            beta1 = (13.0_wp/12.0_wp)*(v1 - 2.0_wp*v2 + v3)**2 &
                  + (1.0_wp/4.0_wp)*(v1 - v3)**2
            beta2 = (13.0_wp/12.0_wp)*(v2 - 2.0_wp*v3 + v4)**2 &
                  + (1.0_wp/4.0_wp)*(3.0_wp*v2 - 4.0_wp*v3 + v4)**2
            
            alpha0 = 0.1_wp / (this%epsilon + beta0)**2
            alpha1 = 0.6_wp / (this%epsilon + beta1)**2
            alpha2 = 0.3_wp / (this%epsilon + beta2)**2
            alpha = alpha0 + alpha1 + alpha2
            w0 = alpha0 / alpha
            w1 = alpha1 / alpha
            w2 = alpha2 / alpha
            
            q0 = (1.0_wp/3.0_wp)*v0 - (7.0_wp/6.0_wp)*v1 + (11.0_wp/6.0_wp)*v2
            q1 = (-1.0_wp/6.0_wp)*v1 + (5.0_wp/6.0_wp)*v2 + (1.0_wp/3.0_wp)*v3
            q2 = (1.0_wp/3.0_wp)*v2 + (5.0_wp/6.0_wp)*v3 - (1.0_wp/6.0_wp)*v4
            
            qL(i-2) = w0 * q0 + w1 * q1 + w2 * q2
            
            ! 计算右界面值 qR(i-2)  
            ! （使用对称模板）
            beta0 = (13.0_wp/12.0_wp)*(v0 - 2.0_wp*v1 + v2)**2 &
                  + (1.0_wp/4.0_wp)*(v0 - 4.0_wp*v1 + 3.0_wp*v2)**2
            beta1 = (13.0_wp/12.0_wp)*(v1 - 2.0_wp*v2 + v3)**2 &
                  + (1.0_wp/4.0_wp)*(v1 - v3)**2
            beta2 = (13.0_wp/12.0_wp)*(v2 - 2.0_wp*v3 + v4)**2 &
                  + (1.0_wp/4.0_wp)*(3.0_wp*v2 - 4.0_wp*v3 + v4)**2
            
            alpha0 = 0.3_wp / (this%epsilon + beta0)**2
            alpha1 = 0.6_wp / (this%epsilon + beta1)**2
            alpha2 = 0.1_wp / (this%epsilon + beta2)**2
            alpha = alpha0 + alpha1 + alpha2
            w0 = alpha0 / alpha
            w1 = alpha1 / alpha
            w2 = alpha2 / alpha
            
            q0 = (-1.0_wp/6.0_wp)*v0 + (5.0_wp/6.0_wp)*v1 + (1.0_wp/3.0_wp)*v2
            q1 = (1.0_wp/3.0_wp)*v1 + (5.0_wp/6.0_wp)*v2 - (1.0_wp/6.0_wp)*v3
            q2 = (11.0_wp/6.0_wp)*v2 - (7.0_wp/6.0_wp)*v3 + (1.0_wp/3.0_wp)*v4
            
            qR(i-2) = w0 * q0 + w1 * q1 + w2 * q2
        end do
    end subroutine weno5_reconstruct
    
    subroutine weno5_info(this)
        class(weno5_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  类型: WENO-5重构器"
    end subroutine weno5_info
    
end module weno5_reconstructor_module