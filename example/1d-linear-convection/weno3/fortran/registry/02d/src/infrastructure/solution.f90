! src/infrastructure/solution.f90
module solution_module
    use base_modules, only: wp, ip
    use domain_module, only: domain_type
    
    implicit none
    private
    public :: wp, ip, solution_type, solution_create, solution_reset
    
    type :: solution_type
        type(domain_type), pointer :: domain => null()
        real(wp), allocatable :: u(:)           ! 当前解（包含ghost）
        real(wp), allocatable :: un(:)          ! 旧解
        real(wp), allocatable :: q_face_left(:) ! 左界面值
        real(wp), allocatable :: q_face_right(:)! 右界面值
        real(wp), allocatable :: flux(:)        ! 通量
        real(wp), allocatable :: res(:)         ! 残差
    contains
        procedure :: initialize => solution_initialize
        procedure :: update_old_field => solution_update_old_field
        procedure :: print_info => solution_print_info
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
        
        allocate(solution%u(ntcells), source=0.0_wp)
        allocate(solution%un(ntcells), source=0.0_wp)
        allocate(solution%q_face_left(nnodes), source=0.0_wp)
        allocate(solution%q_face_right(nnodes), source=0.0_wp)
        allocate(solution%flux(nnodes), source=0.0_wp)
        allocate(solution%res(ncells), source=0.0_wp)
        
        if (domain%config%verbose) then
            print *, "[SOLUTION] Created with arrays:"
            print *, "           u size: ", size(solution%u)
            print *, "           flux size: ", size(solution%flux)
            print *, "           res size: ", size(solution%res)
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
        
        ! 同步旧场
        call this%update_old_field()
        
        if (domain%config%verbose) then
            print *, "[SOLUTION] Initialized"
        end if
    end subroutine solution_initialize
    
    subroutine solution_update_old_field(this)
        class(solution_type), intent(inout) :: this
        this%un = this%u
    end subroutine solution_update_old_field
    
    subroutine solution_reset(this)
        class(solution_type), intent(inout) :: this
        this%u = 0.0_wp
        this%un = 0.0_wp
        this%q_face_left = 0.0_wp
        this%q_face_right = 0.0_wp
        this%flux = 0.0_wp
        this%res = 0.0_wp
    end subroutine solution_reset
    
    subroutine solution_print_info(this)
        class(solution_type), intent(in) :: this
        
        print *, "=== Solution Information ==="
        print *, "u size: ", size(this%u)
        print *, "un size: ", size(this%un)
        print *, "q_face_left size: ", size(this%q_face_left)
        print *, "q_face_right size: ", size(this%q_face_right)
        print *, "flux size: ", size(this%flux)
        print *, "res size: ", size(this%res)
        if (allocated(this%u)) then
            print *, "u range: ", minval(this%u), " to ", maxval(this%u)
        end if
        print *, "============================"
    end subroutine solution_print_info

end module solution_module