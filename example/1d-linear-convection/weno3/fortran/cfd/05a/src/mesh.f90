! ==================== mesh.f90 ====================
! Mesh module for OneFLOW-CFD

module mesh_module
    use kinds, only: dp
    implicit none
    
    private
    
    ! ===================================================================
    ! Mesh Type
    ! ===================================================================
    type, public :: MeshType
        real(dp) :: xmin = 0.0_dp
        real(dp) :: xmax = 2.0_dp
        !integer :: ncells = 40
		integer :: ncells = 200
        integer :: nnodes = 0
        integer :: nx = 0
        real(dp) :: L = 0.0_dp
        real(dp) :: dx = 0.0_dp
        real(dp), allocatable :: x(:), xcc(:)
        logical :: initialized = .false.
    contains
        procedure :: init => mesh_init
        procedure :: print_info => mesh_print_info
    end type MeshType
    
contains
    
    ! ===================================================================
    ! Mesh Initialization
    ! ===================================================================
    
 ! ===================================================================
    ! Mesh Initialization
    ! ===================================================================
    subroutine mesh_init(this, xmin, xmax, ncells)
        class(MeshType), intent(inout) :: this
        real(dp), intent(in), optional :: xmin, xmax
        integer, intent(in), optional :: ncells
        
        integer :: i
        
        ! 1. 更新参数（如果提供）
        if (present(xmin)) this%xmin = xmin
        if (present(xmax)) this%xmax = xmax
        if (present(ncells)) this%ncells = ncells
        
        ! 2. 验证参数
        if (this%xmax <= this%xmin) then
            error stop "Error: xmax must be greater than xmin"
        end if
        
        if (this%ncells <= 0) then
            error stop "Error: ncells must be positive"
        end if
        
        ! 3. 计算派生参数
        this%nnodes = this%ncells + 1
        this%L = this%xmax - this%xmin
        this%dx = this%L / real(this%ncells, dp)
        
        ! 4. 分配内存
        if (allocated(this%x)) deallocate(this%x)
        if (allocated(this%xcc)) deallocate(this%xcc)
        
        allocate(this%x(this%nnodes))
        allocate(this%xcc(this%ncells))
        
        ! 5. 生成网格
        ! 节点坐标
        do i = 1, this%nnodes
            this%x(i) = this%xmin + (i-1) * this%dx
        end do
        
        ! 单元中心坐标
        do i = 1, this%ncells
            this%xcc(i) = 0.5_dp * (this%x(i) + this%x(i+1))
        end do
        
        this%initialized = .true.
        
        ! 6. 打印信息
        call this%print_info()
    end subroutine mesh_init
    
    subroutine mesh_print_info(this)
        class(MeshType), intent(in) :: this
        
        if (.not. this%initialized) then
            print *, "Mesh: Not initialized"
            return
        end if
        
        print *, "=== Mesh Information ==="
        print *, "Domain:    [", this%xmin, ", ", this%xmax, "]"
        print *, "Length:     ", this%L
        print *, "Cells:      ", this%ncells
        print *, "Nodes:      ", this%nnodes
        print *, "Δx:         ", this%dx
        print *, "Cell aspect: ", this%dx / this%L
        print *, ""
    end subroutine    
    
end module mesh_module