! ==================== weno3.f90 ================================
! WENO (Weighted Essentially Non-Oscillatory) reconstructor for OneFLOW-CFD
!
! Purpose: Implement WENO reconstruction scheme
! Reference: Jiang & Shu (1996)
! Author: OneFLOW-CFD Team
! Date: Created
! ===============================================================

module weno3_module
    use kinds, only: dp
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    implicit none
    
    private
    
    ! WENO重构器类型
    type, public :: WenoReconstructorType
    contains
        procedure :: reconstruct => weno_reconstruct
    end type WenoReconstructorType
    
    ! WENO相关参数
    real(dp), parameter, public :: eps_weno = 1.0e-6_dp
    
contains
    
    ! WENO reconstruction
    subroutine weno_reconstruct(this, q, domain, solution)
        class(WenoReconstructorType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        call weno3L_periodic(q, solution%q_face_left, domain, solution)
        call weno3R_periodic(q, solution%q_face_right, domain, solution)
    end subroutine weno_reconstruct
	
    ! WENO-3 reconstruction for left interface
    subroutine weno3L_periodic(u, qL, domain, solution)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: qL(:)
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        integer :: i, j, nghosts, ist, ied
        
        nghosts = domain%nghosts
        ist = domain%ist
        ied = domain%ied
        
        do i = ist-1, ied-1
            j = i - (ist - 1) + 1
            qL(j) = wc3L(u(i-1), u(i), u(i+1))
        end do
        
    end subroutine weno3L_periodic
    
    ! WENO-3 reconstruction for right interface
    subroutine weno3R_periodic(u, qR, domain, solution)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: qR(:)
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        integer :: i, j, nghosts, ist, ied
        
        nghosts = domain%nghosts
        ist = domain%ist
        ied = domain%ied
        
        do i = ist, ied
            j = i - ist + 1
            qR(j) = wc3R(u(i-1), u(i), u(i+1))
        end do
        
    end subroutine weno3R_periodic
    
    ! WENO-3 nonlinear weights for left-biased stencil
    function wc3L(v1, v2, v3) result(f)
        real(dp), intent(in) :: v1, v2, v3
        real(dp) :: f, s0, s1, d0, d1, c0, c1, w0, w1, q0, q1
        
        ! Smoothness indicators
        s0 = (v3 - v2)**2
        s1 = (v2 - v1)**2
        
        ! Nonlinear weights
        d0 = 2.0_dp/3.0_dp
        d1 = 1.0_dp/3.0_dp
        
        c0 = d0 / ((eps_weno + s0)**2)
        c1 = d1 / ((eps_weno + s1)**2)
        
        w0 = c0 / (c0 + c1)
        w1 = c1 / (c0 + c1)
        
        ! Candidate stencils
        q0 = 0.5_dp * v2 + 0.5_dp * v3
        q1 = -0.5_dp * v1 + 1.5_dp * v2
        
        ! Reconstructed value
        f = w0 * q0 + w1 * q1
    end function wc3L
    
    ! WENO-3 nonlinear weights for right-biased stencil
    function wc3R(v1, v2, v3) result(f)
        real(dp), intent(in) :: v1, v2, v3
        real(dp) :: f, s0, s1, d0, d1, c0, c1, w0, w1, q0, q1
        
        ! Smoothness indicators
        s0 = (v2 - v1)**2
        s1 = (v3 - v2)**2
        
        ! Nonlinear weights
        d0 = 2.0_dp/3.0_dp
        d1 = 1.0_dp/3.0_dp
        
        c0 = d0 / ((eps_weno + s0)**2)
        c1 = d1 / ((eps_weno + s1)**2)
        
        w0 = c0 / (c0 + c1)
        w1 = c1 / (c0 + c1)
        
        ! Candidate stencils
        q0 = 0.5_dp * v1 + 0.5_dp * v2
        q1 = 1.5_dp * v2 - 0.5_dp * v3
        
        ! Reconstructed value
        f = w0 * q0 + w1 * q1
    end function wc3R	
    
end module weno3_module