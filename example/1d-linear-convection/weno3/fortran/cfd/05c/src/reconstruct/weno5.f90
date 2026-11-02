! ==================== weno5.f90 ================================
! WENO5 (5th order Weighted Essentially Non-Oscillatory) reconstructor for OneFLOW-CFD
!
! Purpose: Implement WENO5 reconstruction scheme
! Reference: Jiang & Shu (1996), Shu (1998)
! Author: OneFLOW-CFD Team
! Date: Created
! ===============================================================

module weno5_module
    use kinds, only: dp
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    implicit none
    
    private
    
    ! WENO5重构器类型
    type, public :: Weno5ReconstructorType
    contains
        procedure :: reconstruct => weno5_reconstruct
    end type Weno5ReconstructorType
    
    ! WENO相关参数
    real(dp), parameter, public :: eps = 1.0e-6_dp
    
contains
    
    ! WENO5 reconstruction
    subroutine weno5_reconstruct(this, q, domain, solution)
        class(Weno5ReconstructorType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        call weno5L(q, solution%q_face_left, domain)
        call weno5R(q, solution%q_face_right, domain)
    end subroutine weno5_reconstruct
    
    ! WENO5 reconstruction for left interface
    subroutine weno5L(u, qL, domain)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: qL(:)
        type(ComputationalDomainType), intent(in) :: domain
        
        integer :: i, j, nghosts, ist, ied
        
        nghosts = domain%nghosts
        ist = domain%ist
        ied = domain%ied
        
        ! 周期边界处理
        do i = ist-1, ied-1
            j = i - (ist - 1) + 1
            
            ! 获取模板数据（考虑周期边界）
            qL(j) = weno5L_point(u(i-2),u(i-1),u(i  ),u(i+1),u(i+2))
        end do
        
    end subroutine weno5L
    
    ! WENO5 reconstruction for right interface
    subroutine weno5R(u, qR, domain)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: qR(:)
        type(ComputationalDomainType), intent(in) :: domain
        
        integer :: i, j, nghosts, ist, ied
        
        nghosts = domain%nghosts
        ist = domain%ist
        ied = domain%ied
        
        ! 周期边界处理
        do i = ist, ied
            j = i - ist + 1
            
            ! 获取模板数据（考虑周期边界）
			qR(j) = weno5R_point(u(i-2),u(i-1),u(i  ),u(i+1),u(i+2))
        end do
        
    end subroutine weno5R
    
    ! WENO5 reconstruction at a point for left interface
    function weno5L_point(v1, v2, v3, v4, v5) result(f)
        real(dp), intent(in) :: v1, v2, v3, v4, v5
        real(dp) :: f
        
        real(dp) :: d0, d1, d2
        real(dp) :: s0, s1, s2
		real(dp) :: q0, q1, q2
        real(dp) :: alpha0, alpha1, alpha2, alpha
		real(dp) :: beta0, beta1, beta2
        real(dp) :: w0, w1, w2
       
        ! 光滑度指示器
        beta0 = (13.0_dp/12.0_dp)*(v1 - 2.0_dp*v2 + v3)**2 + (1.0_dp/4.0_dp)*(v1 - 4.0_dp*v2 + 3.0_dp*v3)**2
        beta1 = (13.0_dp/12.0_dp)*(v2 - 2.0_dp*v3 + v4)**2 + (1.0_dp/4.0_dp)*(v2 - v4)**2
        beta2 = (13.0_dp/12.0_dp)*(v3 - 2.0_dp*v4 + v5)**2 + (1.0_dp/4.0_dp)*(3.0_dp*v3 - 4.0_dp*v4 + v5)**2
             
		d0 = 1.0_dp/10.0_dp
		d1 = 3.0_dp/5.0_dp
		d2 = 3.0_dp/10.0_dp
		alpha0 = d0 / (eps + beta0)**2
		alpha1 = d1 / (eps + beta1)**2
		alpha2 = d2 / (eps + beta2)**2
        
		alpha = alpha0 + alpha1 + alpha2
        
        w0 = alpha0 / alpha
        w1 = alpha1 / alpha
        w2 = alpha2 / alpha

		q0 = 1.0_dp/3.0_dp*v1 - 7.0_dp/6.0_dp*v2 + 11.0_dp/6.0_dp*v3  ! r=2
		q1 = -1.0_dp/6.0_dp*v2 + 5.0_dp/6.0_dp*v3 + 1.0_dp/3.0_dp*v4  ! r=1
		q2 = 1.0_dp/3.0_dp*v3 + 5.0_dp/6.0_dp*v4 - 1.0_dp/6.0_dp*v5   ! r=0
        
        ! 加权重构
        f = w0 * q0 + w1 * q1 + w2 * q2
        
    end function weno5L_point
    
    ! WENO5 reconstruction at a point for right interface
    function weno5R_point(v1, v2, v3, v4, v5) result(f)
        real(dp), intent(in) :: v1, v2, v3, v4, v5
        real(dp) :: f
        
        real(dp) :: d0, d1, d2
        real(dp) :: s0, s1, s2
		real(dp) :: q0, q1, q2
        real(dp) :: alpha0, alpha1, alpha2, alpha
		real(dp) :: beta0, beta1, beta2
        real(dp) :: w0, w1, w2
        
		beta0 = (13.0_dp/12.0_dp)*(v1 - 2*v2 + v3)**2 + (1.0_dp/4.0_dp)*(v1 - 4*v2 + 3*v3)**2
		beta1 = (13.0_dp/12.0_dp)*(v2 - 2*v3 + v4)**2 + (1.0_dp/4.0_dp)*(v2 - v4)**2
		beta2 = (13.0_dp/12.0_dp)*(v3 - 2*v4 + v5)**2 + (1.0_dp/4.0_dp)*(3*v3 - 4*v4 + v5)**2
		
		d0 = 3.0_dp/10.0_dp
		d1 = 3.0_dp/5.0_dp
		d2 = 1.0_dp/10.0_dp
		alpha0 = d0 / (eps + beta0)**2
		alpha1 = d1 / (eps + beta1)**2
		alpha2 = d2 / (eps + beta2)**2
		alpha = alpha0 + alpha1 + alpha2
		w0 = alpha0 / alpha
		w1 = alpha1 / alpha
		w2 = alpha2 / alpha
		
		q0 = -1.0_dp/6.0_dp*v1 + 5.0_dp/6.0_dp*v2 + 1.0_dp/3.0_dp*v3  ! r=2
		q1 = 1.0_dp/3.0_dp*v2 + 5.0_dp/6.0_dp*v3 - 1.0_dp/6.0_dp*v4  ! r=1
		q2 = 11.0_dp/6.0_dp*v3 - 7.0_dp/6.0_dp*v4 + 1.0_dp/3.0_dp*v5  ! r=0
        
        ! 加权重构
        f = w0 * q0 + w1 * q1 + w2 * q2
        
    end function weno5R_point
    
end module weno5_module