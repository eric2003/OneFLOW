! ==================== reconstruct.f90 ====================
! Configuration module for OneFLOW-CFD

module recon_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    implicit none
    
    private
	public :: init_coef, eno_reconstruct, weno_reconstruct
	
    ! ===================================================================
    ! Solution Type
    ! ===================================================================
    type, public :: ReconType
        type(EnoReconstructorType), allocatable :: eno_recon
        type(WenoReconstructorType), allocatable :: weno_recon
    contains
        !procedure :: init => recon_init
    end type ReconType
	
    ! ENO重构器：独立类型（去掉 extends(ReconstructorType)）
    type, public :: EnoReconstructorType
        integer :: spatial_order
        integer :: ntcells
        integer, allocatable :: lmc(:)
        real(dp), allocatable :: coef(:,:)
        real(dp), allocatable :: dd(:,:)
    contains
        procedure :: reconstruct => eno_reconstruct
        procedure :: init => eno_init
    end type EnoReconstructorType
    
    ! WENO重构器：独立类型（去掉 extends(ReconstructorType)）
    type, public :: WenoReconstructorType
    contains
        procedure :: reconstruct => weno_reconstruct
    end type WenoReconstructorType
    
    ! ===================================================================
    ! Module Variables
    ! ===================================================================
    real(dp), parameter :: eps_weno = 1.0e-6_dp    
    
contains

	! Initialize ENO/WENO coefficients
	subroutine init_coef(spatial_order, coef)
		!use, intrinsic :: iso_fortran_env, only: dp => real64
		integer, intent(in) :: spatial_order
		real(dp), intent(out) :: coef(:,:)
		
		coef = 0.0_dp
		
		select case(spatial_order)
		case(1)
			coef(1,1) = 1.0_dp
			coef(2,1) = 1.0_dp
			
		case(2)
			coef(1,1:2) = [ 3.0_dp/2.0_dp, -1.0_dp/2.0_dp ]
			coef(2,1:2) = [ 1.0_dp/2.0_dp,  1.0_dp/2.0_dp ]
			coef(3,1:2) = [ -1.0_dp/2.0_dp, 3.0_dp/2.0_dp ]
			
		case(3)
			coef(1,1:3) = [ 11.0_dp/6.0_dp, -7.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
			coef(2,1:3) = [  1.0_dp/3.0_dp,  5.0_dp/6.0_dp, -1.0_dp/6.0_dp ]
			coef(3,1:3) = [ -1.0_dp/6.0_dp,  5.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
			coef(4,1:3) = [  1.0_dp/3.0_dp, -7.0_dp/6.0_dp, 11.0_dp/6.0_dp ]
			
		case default
			error stop "Unsupported spatial order"
		end select
    end subroutine init_coef
    
    ! ENO reconstructor initialization
    subroutine eno_init(this, spatial_order, ntcells)
        class(EnoReconstructorType), intent(inout) :: this
        integer, intent(in) :: spatial_order
        integer, intent(in) :: ntcells
        
        this%spatial_order = spatial_order
        this%ntcells = ntcells
        
        allocate(this%lmc(ntcells))
        allocate(this%coef(spatial_order+1, spatial_order))
        allocate(this%dd(spatial_order, ntcells))
        
        this%lmc = 0
        this%coef = 0.0_dp
        this%dd = 0.0_dp
        
        ! Initialize coefficients
		call init_coef(spatial_order, this%coef)

		print *, "ENO reconstructor initialized:"
		print *, "  spatial_order = ", spatial_order
		print *, "  ntcells = ", ntcells
    end subroutine eno_init
    
   ! ENO reconstruction
    subroutine eno_reconstruct(this, q, domain, solution)
        class(EnoReconstructorType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        integer :: i, j, m, k1, k2, r1, r2
        integer :: nghosts, ist, ied
		
        nghosts = domain%nghosts
        ist = domain%ist
        ied = domain%ied

		! 1. 差商计算 (dd[1,:] = q)
		this%dd(1, :) = q(:)
		do m = 2, this%spatial_order
			do j = 1, domain%ntcells - m + 1
				this%dd(m, j) = this%dd(m-1, j+1) - this%dd(m-1, j)
			end do
		end do
		
		! 2. 选择 smoothest stencil
		do i = ist - 1, ied
			this%lmc(i) = i
			do m = 2, this%spatial_order
				if ( abs(this%dd(m, this%lmc(i) - 1) ) < abs(this%dd(m, this%lmc(i)))) then
					this%lmc(i) = this%lmc(i) - 1
                end if
            end do
        end do
		
		! 3. 重构界面值
		do i = ist, ied
			j = i - ist + 1  ! 1-based index for interfaces
			k1 = this%lmc(i - 1)
			k2 = this%lmc(i)
			r1 = (i - 1) - k1 + 1
			r2 = i - k2 + 1
			solution%q_face_left(j) = 0.0_dp
			solution%q_face_right(j) = 0.0_dp
			do m = 1, this%spatial_order
				solution%q_face_left(j) = solution%q_face_left(j) + q(k1 + m - 1) * this%coef(r1 + 1, m)
				solution%q_face_right(j) = solution%q_face_right(j) + q(k2 + m - 1) * this%coef(r2, m)
			end do
		end do
       
    end subroutine eno_reconstruct  
    
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

end module recon_module