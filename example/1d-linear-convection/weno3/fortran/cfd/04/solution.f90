! ==================== solution.f90 ====================
! Solution variables module for OneFLOW-CFD

module solution_module
    use kinds, only: dp
    use domain_module, only: ComputationalDomainType
    implicit none
    
    private
	public :: update_oldfield
    
    ! ===================================================================
    ! Solution Type
    ! ===================================================================
    type, public :: SolutionType
        type(ComputationalDomainType) :: domain
        real(dp), allocatable :: q_face_left(:), q_face_right(:)
        real(dp), allocatable :: flux(:), res(:)
        real(dp), allocatable :: u(:), un(:)
    contains
        procedure :: init => solution_init
    end type SolutionType
    
    ! 类型已经在声明时标记为public，不需要额外声明
    
contains
    
    ! ===================================================================
    ! Solution Initialization
    ! ===================================================================
    
    ! Solution initialization
    subroutine solution_init(this, domain)
        class(SolutionType), intent(inout) :: this
        type(ComputationalDomainType), intent(in) :: domain
        
        this%domain = domain
        
        allocate(this%q_face_left(domain%mesh%nnodes))
        allocate(this%q_face_right(domain%mesh%nnodes))
        allocate(this%flux(domain%mesh%nnodes))
        allocate(this%res(domain%mesh%ncells))
        allocate(this%u(domain%ntcells))
        allocate(this%un(domain%ntcells))
        
        this%q_face_left = 0.0_dp
        this%q_face_right = 0.0_dp
        this%flux = 0.0_dp
        this%res = 0.0_dp
        this%u = 0.0_dp
        this%un = 0.0_dp
    end subroutine solution_init
	
    ! Update old field
    subroutine update_oldfield(qn, q, n)
        real(dp), intent(out) :: qn(:)
        real(dp), intent(in) :: q(:)
        integer, intent(in) :: n
        
        qn(1:n) = q(1:n)
    end subroutine update_oldfield	
    
end module solution_module