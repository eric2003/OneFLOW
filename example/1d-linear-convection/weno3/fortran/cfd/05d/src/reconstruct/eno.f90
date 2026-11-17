! ==================== eno.f90 =================================
! ENO (Essentially Non-Oscillatory) reconstructor for OneFLOW-CFD
!
! Purpose: Implement ENO reconstruction scheme
! Reference: Harten, Engquist, Osher, Chakravarthy (1987)
! Author: OneFLOW-CFD Team
! Date: Created
! ==============================================================

module eno_module
    use kinds, only: dp
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    use recon_coefficients_module, only: init_coef
    implicit none
    
    private
   
    ! ENO重构器类型
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
    
contains
    
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
    
end module eno_module