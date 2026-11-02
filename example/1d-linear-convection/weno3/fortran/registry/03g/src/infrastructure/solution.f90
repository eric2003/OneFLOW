! src/infrastructure/solution.f90（修正版）
module solution_module
    use base_modules, only: wp, ip
    use domain_module, only: domain_type
    
    implicit none
    private
    public :: wp, ip, solution_type, solution_create, solution_reset
    
    type :: solution_type
        type(domain_type), pointer :: domain => null()
        real(wp), allocatable :: u(:)           ! 当前解（含ghost）
        real(wp), allocatable :: un(:)          ! 旧解
        real(wp), allocatable :: q_face_left(:) ! 左界面值
        real(wp), allocatable :: q_face_right(:)! 右界面值
        real(wp), allocatable :: flux(:)        ! 通量
        real(wp), allocatable :: res(:)         ! 残差
    contains
        procedure :: initialize => solution_initialize
        procedure :: update_old_field => solution_update_old_field
        procedure :: print_info => solution_print_info
        procedure :: reset => solution_reset_instance
    end type solution_type
    
contains

    function solution_create(domain) result(solution)
        type(domain_type), target, intent(in) :: domain
        type(solution_type) :: solution
        
        integer(ip) :: ncells, nnodes, ntcells
        
        solution%domain => domain
        
        ncells = domain%mesh%ncells
        nnodes = domain%mesh%nnodes
        ntcells = domain%ntcells
        
        ! 分配数组（与Julia solution.jl一致）
        allocate(solution%u(ntcells), source=0.0_wp)
        allocate(solution%un(ntcells), source=0.0_wp)
        allocate(solution%q_face_left(nnodes), source=0.0_wp)
        allocate(solution%q_face_right(nnodes), source=0.0_wp)
        allocate(solution%flux(nnodes), source=0.0_wp)
        allocate(solution%res(ncells), source=0.0_wp)
        
        if (domain%config%verbose) then
            print *, "[SOLUTION] Created:"
            print *, "  u size: ", size(solution%u), " (with ghosts)"
            print *, "  flux size: ", size(solution%flux)
            print *, "  res size: ", size(solution%res)
        end if
    end function solution_create
    
    subroutine solution_initialize(this, initial_values)
        class(solution_type), intent(inout) :: this
        real(wp), intent(in), optional :: initial_values(:)
        
        integer(ip) :: i, idx
        type(domain_type), pointer :: domain
        
        domain => this%domain
        
        if (present(initial_values)) then
            ! 应用初始值到物理区域
            do i = domain%ist, domain%ied - 1
                idx = i - domain%ist + 1
                if (idx <= size(initial_values)) then
                    this%u(i) = initial_values(idx)
                end if
            end do
        else
            ! 默认为0
            this%u = 0.0_wp
        end if
        
        ! 同步旧场（与Julia的update_old_field一致）
        call this%update_old_field()
        
        if (domain%config%verbose) then
            print *, "[SOLUTION] Initialized"
            print *, "  u range: ", minval(this%u), " to ", maxval(this%u)
        end if
    end subroutine solution_initialize
    
    subroutine solution_update_old_field(this)
        class(solution_type), intent(inout) :: this
        this%un = this%u  ! 与Julia的 un .= u 一致
    end subroutine solution_update_old_field
    
    subroutine solution_reset_instance(this)
        class(solution_type), intent(inout) :: this
        call solution_reset(this)
    end subroutine solution_reset_instance
    
    subroutine solution_reset(solution)
        type(solution_type), intent(inout) :: solution
        
        if (allocated(solution%u)) solution%u = 0.0_wp
        if (allocated(solution%un)) solution%un = 0.0_wp
        if (allocated(solution%q_face_left)) solution%q_face_left = 0.0_wp
        if (allocated(solution%q_face_right)) solution%q_face_right = 0.0_wp
        if (allocated(solution%flux)) solution%flux = 0.0_wp
        if (allocated(solution%res)) solution%res = 0.0_wp
        
        if (associated(solution%domain) .and. solution%domain%config%verbose) then
            print *, "[SOLUTION] Reset"
        end if
    end subroutine solution_reset
    
    subroutine solution_print_info(this)
        class(solution_type), intent(in) :: this
        
        print *, "=== Solution Information ==="
        print *, "Arrays:"
        print *, "  u: ", size(this%u), " elements"
        print *, "  un: ", size(this%un), " elements"
        print *, "  q_face_left: ", size(this%q_face_left), " elements"
        print *, "  q_face_right: ", size(this%q_face_right), " elements"
        print *, "  flux: ", size(this%flux), " elements"
        print *, "  res: ", size(this%res), " elements"
        
        if (allocated(this%u)) then
            print *, "Values:"
            print *, "  u min/max: ", minval(this%u), maxval(this%u)
            print *, "  un min/max: ", minval(this%un), maxval(this%un)
        end if
        print *, "============================"
    end subroutine solution_print_info

end module solution_module