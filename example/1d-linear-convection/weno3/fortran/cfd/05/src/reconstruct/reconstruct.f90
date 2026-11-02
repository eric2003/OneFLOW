! ==================== reconstruct.f90 ====================
! Main reconstruction module for OneFLOW-CFD
!
! Purpose: Unified interface for reconstruction schemes
! Author: OneFLOW-CFD Team
! Date: Created
! ========================================================

module recon_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    
    ! 导入具体重构器模块
    use eno_module, only: EnoReconstructorType
    use weno3_module, only: WenoReconstructorType
    
    implicit none
    
    private
    public :: reconstruction
    
    ! 重构器组合类型
    type, public :: ReconType
        type(EnoReconstructorType), allocatable :: eno_recon
        type(WenoReconstructorType), allocatable :: weno_recon
    contains
        procedure :: init => recon_init
        procedure :: reconstruct => recon_reconstruct
    end type ReconType
    
contains
    
    ! 重构器初始化
    subroutine recon_init(this, config, domain)
        class(ReconType), intent(inout) :: this
        type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        
        if (trim(config%recon_scheme) == "eno") then
            ! 分配并初始化ENO重构器
            allocate(this%eno_recon)
            call this%eno_recon%init(config%spatial_order, domain%ntcells)
        else if (trim(config%recon_scheme) == "weno") then
            ! 分配WENO重构器
            allocate(this%weno_recon)
        else
            error stop "Unknown reconstruction scheme: "//trim(config%recon_scheme)
        end if
    end subroutine recon_init
    
    ! 重构器执行
    subroutine recon_reconstruct(this, q, config, domain, solution)
        class(ReconType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        character(len=:), allocatable :: trimmed_scheme
        
        trimmed_scheme = trim(config%recon_scheme)
        
        if (trimmed_scheme == "eno") then
            call this%eno_recon%reconstruct(q, domain, solution)
        else if (trimmed_scheme == "weno") then
            call this%weno_recon%reconstruct(q, domain, solution)
        else
            error stop "Unknown reconstruction scheme in recon_reconstruct(): "//trimmed_scheme
        end if
        
    end subroutine recon_reconstruct
    
    ! 重构函数（兼容性接口）
    subroutine reconstruction(this, q, config, domain, solution)
        type(ReconType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        call this%reconstruct(q, config, domain, solution)
    end subroutine reconstruction
    
end module recon_module
