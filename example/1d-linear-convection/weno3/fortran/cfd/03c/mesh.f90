! ==================== mesh.f90 ====================
! Mesh module for OneFLOW-CFD

module mesh_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    private
    
    ! ===================================================================
    ! Mesh Type
    ! ===================================================================
    type, public :: MeshType
        real(dp) :: xmin = 0.0_dp
        real(dp) :: xmax = 2.0_dp
        integer :: ncells = 40
        integer :: nnodes = 0
        integer :: nx = 0
        real(dp) :: L = 0.0_dp
        real(dp) :: dx = 0.0_dp
        real(dp), allocatable :: x(:), xcc(:)
    contains
        procedure :: init => mesh_init
    end type MeshType
    
contains
    
    ! ===================================================================
    ! Mesh Initialization
    ! ===================================================================
    
    ! Mesh initialization
    subroutine mesh_init(this)
        class(MeshType), intent(inout) :: this
        integer :: i
        
        this%nnodes = this%ncells + 1
        this%nx = this%ncells
        this%L = this%xmax - this%xmin
        this%dx = this%L / real(this%ncells, dp)
        
        allocate(this%x(this%nnodes), this%xcc(this%ncells))
        
        ! Node coordinates
        do i = 1, this%nnodes
            this%x(i) = this%xmin + (i-1) * this%dx
        end do
        
        ! Cell center coordinates
        do i = 1, this%ncells
            this%xcc(i) = 0.5_dp * (this%x(i) + this%x(i+1))
        end do
    end subroutine mesh_init
    
end module mesh_module