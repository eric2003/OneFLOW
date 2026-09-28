! ==================== domain.f90 ====================
! Computational domain module for OneFLOW-CFD

module domain_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use mesh_module, only: MeshType
    implicit none
    
    private
    
    ! ===================================================================
    ! Computational Domain Type
    ! ===================================================================
    type, public :: ComputationalDomainType
        type(MeshType) :: mesh
        type(CfdConfigType), pointer :: config => null()
        integer :: nghosts = 0
        integer :: ist = 0
        integer :: ied = 0
        integer :: ntcells = 0
    contains
        procedure :: init => domain_init
    end type ComputationalDomainType
    
    ! 类型已经在声明时标记为public，不需要额外声明
    
contains
    
    ! ===================================================================
    ! Domain Initialization
    ! ===================================================================
    
    ! Domain initialization
    subroutine domain_init(this, config, mesh)
        class(ComputationalDomainType), intent(inout) :: this
        type(MeshType), intent(in) :: mesh
        type(CfdConfigType), target, intent(in) :: config
        
        this%mesh = mesh
        this%config => config
        
        ! Calculate ghost cells
        if (trim(config%recon_scheme) == "eno") then
            this%nghosts = config%spatial_order
        else if (trim(config%recon_scheme) == "weno") then
            this%nghosts = config%spatial_order / 2 + 1
        else
            error stop "Unknown reconstruction scheme"
        end if
        
        this%ist = this%nghosts + 1
        !this%ied = this%ist + mesh%ncells - 1
        this%ied = this%ist + mesh%ncells
        this%ntcells = mesh%ncells + 2 * this%nghosts
        
        print *, "Domain initialized:"
        print *, "  mesh.ncells = ", mesh%ncells
        print *, "  spatial_order = ", config%spatial_order
        print *, "  nghosts = ", this%nghosts
        print *, "  ist = ", this%ist, ", ied = ", this%ied
        print *, "  dx = ", mesh%dx
    end subroutine domain_init
    
end module domain_module