! src/numerics/reconstructor/weno3.f90
module weno3_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use reconstructor_base_module, only: reconstructor_base, wp
    use registry_module, only: register_component_with_factory
    implicit none
    
    private
    public :: wp, weno3_reconstructor, create_weno3
    
    type, extends(reconstructor_base) :: weno3_reconstructor
        real(wp) :: epsilon = 1.0e-6_wp
    contains
        procedure :: reconstruct => weno3_reconstruct
        procedure :: info => weno3_info
    end type weno3_reconstructor
    
contains
    
    ! Factory function
    subroutine create_weno3(instance)
        class(*), allocatable, intent(out) :: instance
        type(weno3_reconstructor), allocatable :: weno3
        
        allocate(weno3)
        weno3%name = "WENO3"
        weno3%order = 3
        weno3%epsilon = 1.0e-6_wp
        
        call move_alloc(weno3, instance)
    end subroutine create_weno3
    
    ! WENO-3 reconstruction
    subroutine weno3_reconstruct(this, q, qL, qR)
        class(weno3_reconstructor), intent(in) :: this
        real(wp), intent(in) :: q(:)
        real(wp), intent(out) :: qL(:)
        real(wp), intent(out) :: qR(:)
        
        integer :: i, n
        real(wp) :: beta0, beta1, alpha0, alpha1, alpha, w0, w1
        real(wp) :: q0, q1
        
        n = size(q) - 2  ! Need 2 ghost cells
        
        if (n /= size(qL) .or. n /= size(qR)) then
            print *, "[ERROR] Array size mismatch in weno3_reconstruct"
            return
        end if
        
        do i = 2, n+1
            ! Left interface value qL(i-1)
            beta0 = (q(i) - q(i-1))**2
            beta1 = (q(i+1) - q(i))**2
            
            alpha0 = 1.0_wp/3.0_wp / (this%epsilon + beta0)**2
            alpha1 = 2.0_wp/3.0_wp / (this%epsilon + beta1)**2
            alpha = alpha0 + alpha1
            w0 = alpha0 / alpha
            w1 = alpha1 / alpha
            
            q0 = -0.5_wp*q(i-1) + 1.5_wp*q(i)
            q1 = 0.5_wp*q(i) + 0.5_wp*q(i+1)
            
            qL(i-1) = w0 * q0 + w1 * q1
            
            ! Right interface value qR(i-1)
            beta0 = (q(i) - q(i-1))**2
            beta1 = (q(i+1) - q(i))**2
            
            alpha0 = 2.0_wp/3.0_wp / (this%epsilon + beta0)**2
            alpha1 = 1.0_wp/3.0_wp / (this%epsilon + beta1)**2
            alpha = alpha0 + alpha1
            w0 = alpha0 / alpha
            w1 = alpha1 / alpha
            
            q0 = 0.5_wp*q(i-1) + 0.5_wp*q(i)
            q1 = 1.5_wp*q(i) - 0.5_wp*q(i+1)
            
            qR(i-1) = w0 * q0 + w1 * q1
        end do
    end subroutine weno3_reconstruct
    
    subroutine weno3_info(this)
        class(weno3_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  Type: WENO-3 reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine weno3_info
    
end module weno3_reconstructor_module