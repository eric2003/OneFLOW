! src/infrastructure/mesh.f90
module mesh_module
    use base_modules, only: wp, ip
    
    implicit none
    public :: wp, ip, mesh_type, mesh_init, mesh_print_info
    
    ! 网格类型
    type :: mesh_type
        real(wp) :: xmin = 0.0_wp
        real(wp) :: xmax = 2.0_wp
        integer(ip) :: ncells = 40
        integer(ip) :: nnodes
        integer(ip) :: nx
        real(wp), allocatable :: x(:)      ! 节点坐标
        real(wp), allocatable :: xcc(:)    ! 单元中心坐标
        real(wp) :: L, dx
    contains
        procedure :: init => mesh_init
        procedure :: print_info => mesh_print_info
    end type mesh_type
    
contains

    subroutine mesh_init(this, xmin, xmax, ncells)
        class(mesh_type), intent(inout) :: this
        real(wp), optional, intent(in) :: xmin, xmax
        integer(ip), optional, intent(in) :: ncells
        
        integer(ip) :: i
        
        ! 设置参数
        if (present(xmin)) this%xmin = xmin
        if (present(xmax)) this%xmax = xmax
        if (present(ncells)) this%ncells = ncells
        
        ! 计算
        this%nnodes = this%ncells + 1
        this%nx = this%ncells
        this%L = this%xmax - this%xmin
        this%dx = this%L / real(this%ncells, wp)
        
        ! 分配内存
        if (allocated(this%x)) deallocate(this%x)
        if (allocated(this%xcc)) deallocate(this%xcc)
        
        allocate(this%x(this%nnodes))
        allocate(this%xcc(this%ncells))
        
        ! 生成节点坐标
        do i = 1, this%nnodes
            this%x(i) = this%xmin + (i - 1) * this%dx
        end do
        
        ! 生成单元中心坐标
        do i = 1, this%ncells
            this%xcc(i) = 0.5_wp * (this%x(i) + this%x(i + 1))
        end do
    end subroutine mesh_init
    
    subroutine mesh_print_info(this)
        class(mesh_type), intent(in) :: this
        
        print *, "=== Mesh Information ==="
        print *, "Domain: [", this%xmin, ", ", this%xmax, "]"
        print *, "Cells: ", this%ncells
        print *, "Nodes: ", this%nnodes
        print *, "dx: ", this%dx
        print *, "L: ", this%L
        print *, "========================"
    end subroutine mesh_print_info

end module mesh_module