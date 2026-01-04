! src/infrastructure/mesh.f90
module mesh_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    private
    public :: wp, mesh_type, mesh_init, mesh_print_info

    ! mesh类型
    type :: mesh_type
        real(wp) :: xmin = 0.0_wp
        real(wp) :: xmax = 2.0_wp
        integer :: ncells = 40
        integer :: nnodes
        integer :: nx
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
        integer, optional, intent(in) :: ncells

        integer :: i

        ! Set参数
        if (present(xmin)) this%xmin = xmin
        if (present(xmax)) this%xmax = xmax
        if (present(ncells)) this%ncells = ncells

        ! computation派生参数
        this%nnodes = this%ncells + 1
        this%nx = this%ncells
        this%L = this%xmax - this%xmin
        this%dx = this%L / real(this%ncells, wp)

        ! 分配数组
        if (allocated(this%x)) deallocate(this%x)
        if (allocated(this%xcc)) deallocate(this%xcc)

        allocate(this%x(this%nnodes))
        allocate(this%xcc(this%ncells))

        ! 生成node坐标
        do i = 1, this%nnodes
            this%x(i) = this%xmin + (i - 1) * this%dx
        end do

        ! 生成cell中心坐标
        do i = 1, this%ncells
            this%xcc(i) = 0.5_wp * (this%x(i) + this%x(i + 1))
        end do
    end subroutine mesh_init

    subroutine mesh_print_info(this)
        class(mesh_type), intent(in) :: this

        print *, "=== 网格信息 ==="
        print *, "计算域: [", this%xmin, ", ", this%xmax, "]"
        print *, "单元数: ", this%ncells
        print *, "节点数: ", this%nnodes
        print *, "网格尺寸 dx: ", this%dx
        print *, "域长度 L: ", this%L
        print *, "=========================="
    end subroutine mesh_print_info

end module mesh_module