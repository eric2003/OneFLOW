! src/results.f90 (修正版)
module results_module
    use base_modules, only: wp, ip
    use config_module, only: cfd_config
    use mesh_module, only: mesh_type
    use domain_module, only: domain_type
    use solution_module, only: solution_type
    ! use physics_solver_module, only: physics_solver  ! 暂时注释掉，避免循环依赖
    
    implicit none
    private
    public :: results_saver, results_saver_create, save_results
    
    ! 定义字符串长度常量
    integer, parameter :: STR_LEN = 128
    
    ! 结果类型 - 对应Julia的result字典
    type :: cfd_results
        character(len=STR_LEN) :: solver_name = ""
        real(wp), allocatable :: x(:)           ! 网格坐标（单元中心）
        real(wp), allocatable :: numerical(:)   ! 数值解
        real(wp), allocatable :: analytical(:)  ! 解析解
        character(len=STR_LEN) :: scheme = ""   ! 格式名称
        integer :: order = 0                    ! 阶数
        integer :: rk_order = 0                 ! RK阶数
        real(wp) :: final_time = 0.0_wp         ! 最终时间
        real(wp) :: current_time = 0.0_wp       ! 当前时间
        integer :: total_steps = 0              ! 总步数
        integer :: solver_state = 0             ! 求解器状态
    end type cfd_results
    
    ! 结果保存器
    type :: results_saver
        character(len=STR_LEN) :: base_filename = "results"
        logical :: verbose = .true.
    contains
        procedure :: save_text => results_saver_save_text
        procedure :: save_binary => results_saver_save_binary
        procedure :: load => results_saver_load
    end type results_saver
    
    ! 接口声明
    interface results_saver
        module procedure results_saver_constructor
    end interface
    
contains

    ! 构造函数
    function results_saver_constructor(base_filename, verbose) result(saver)
        character(len=*), optional :: base_filename
        logical, optional :: verbose
        type(results_saver) :: saver
        
        if (present(base_filename)) then
            saver%base_filename = trim(adjustl(base_filename))
        end if
        if (present(verbose)) then
            saver%verbose = verbose
        end if
    end function results_saver_constructor
    
    ! 保持向后兼容的创建函数
    function results_saver_create(base_filename, verbose) result(saver)
        character(len=*), optional :: base_filename
        logical, optional :: verbose
        type(results_saver) :: saver
        
        saver = results_saver_constructor(base_filename, verbose)
    end function results_saver_create
    
    ! 生成文件名（与Julia风格一致）
    function generate_filename(saver, solver_name, mesh_size) result(filename)
        class(results_saver), intent(in) :: saver
        character(len=*), intent(in) :: solver_name
        integer, intent(in) :: mesh_size
        character(len=STR_LEN) :: filename
        
        write(filename, '(A, "_", A, "_", I0, ".dat")') &
            trim(saver%base_filename), trim(solver_name), mesh_size
    end function generate_filename
    
    ! 主保存函数 - 生成与Julia兼容的结果
    subroutine save_results(saver, solver_name, config, mesh, domain, solution, &
                          current_time, total_steps, solver_state)
        class(results_saver), intent(in) :: saver
        character(len=*), intent(in) :: solver_name
        type(cfd_config), intent(in) :: config
        type(mesh_type), intent(in) :: mesh
        type(domain_type), intent(in) :: domain
        type(solution_type), intent(in) :: solution
        real(wp), intent(in) :: current_time
        integer, intent(in) :: total_steps, solver_state
        
        type(cfd_results) :: results
        character(len=STR_LEN) :: filename
        integer :: i, n_physical
        
        ! 准备结果数据
        results%solver_name = trim(solver_name)
        results%scheme = trim(config%recon_scheme)
        results%order = config%spatial_order
        results%rk_order = config%rk_order
        results%final_time = config%final_time
        results%current_time = current_time
        results%total_steps = total_steps
        results%solver_state = solver_state
        
        ! 分配数组
        n_physical = mesh%ncells
        allocate(results%x(n_physical))
        allocate(results%numerical(n_physical))
        allocate(results%analytical(n_physical))
        
        ! 填充网格坐标（单元中心）
        results%x = mesh%xcc
        
        ! 填充数值解（仅物理区域）
        do i = 1, n_physical
            results%numerical(i) = solution%u(domain%ist + i - 1)
        end do
        
        ! 生成解析解（与Julia的exact_solution对应）
        call generate_analytical_solution(results%x, config, results%analytical, current_time)
        
        ! 生成文件名
        filename = generate_filename(saver, solver_name, mesh%ncells)
        
        ! 保存文件
        if (saver%verbose) then
            print *, "[RESULTS] Saving results to: ", trim(filename)
            print *, "  Solver: ", trim(results%solver_name)
            print *, "  Scheme: ", trim(results%scheme), " order ", results%order
            print *, "  Time: ", results%current_time, " / ", results%final_time
            print *, "  Steps: ", results%total_steps
        end if
        
        call saver%save_text(results, filename)
        
        ! 清理
        deallocate(results%x, results%numerical, results%analytical)
    end subroutine save_results
    
    ! 生成解析解（匹配Julia的exact_solution逻辑）
    subroutine generate_analytical_solution(x, config, analytical, current_time)
        real(wp), intent(in) :: x(:), current_time
        type(cfd_config), intent(in) :: config
        real(wp), intent(out) :: analytical(:)
        
        integer :: i, n
        real(wp) :: x_shifted, L
        
        n = size(x)
        L = config%domain_length
        
        select case (trim(config%ic_type))
        case ("step")
            ! 阶跃函数的精确解（周期性）
            do i = 1, n
                ! 周期性平移
                x_shifted = x(i) - config%wave_speed * current_time
                x_shifted = modulo(x_shifted, L)
                if (x_shifted < 0.0_wp) x_shifted = x_shifted + L
                
                ! 阶跃在 [0.5, 1.0] 内为 2.0，其他为 1.0
                if (x_shifted >= 0.5_wp .and. x_shifted <= 1.0_wp) then
                    analytical(i) = 2.0_wp
                else
                    analytical(i) = 1.0_wp
                end if
            end do
            
        case ("sin", "sine")
            ! 正弦波的精确解
            do i = 1, n
                x_shifted = x(i) - config%wave_speed * current_time
                x_shifted = modulo(x_shifted, L)
                analytical(i) = sin(2.0_wp * 3.141592653589793_wp * x_shifted / L)
            end do
            
        case ("gaussian")
            ! 高斯脉冲的精确解
            do i = 1, n
                x_shifted = x(i) - config%wave_speed * current_time
                x_shifted = modulo(x_shifted, L)
                analytical(i) = exp(-50.0_wp * (x_shifted - 1.0_wp)**2)
            end do
            
        case default
            ! 默认：阶跃函数
            do i = 1, n
                x_shifted = x(i) - config%wave_speed * current_time
                x_shifted = modulo(x_shifted, L)
                if (x_shifted >= 0.5_wp .and. x_shifted <= 1.0_wp) then
                    analytical(i) = 2.0_wp
                else
                    analytical(i) = 1.0_wp
                end if
            end do
        end select
    end subroutine generate_analytical_solution
    
    ! 文本格式保存（与Julia的纯文本输出兼容）
    subroutine results_saver_save_text(this, results, filename)
        class(results_saver), intent(in) :: this
        type(cfd_results), intent(in) :: results
        character(len=*), intent(in) :: filename
        
        integer :: i, n, unit, ierr
        
        n = size(results%x)
        
        ! 打开文件
        open(newunit=unit, file=trim(filename), status='replace', &
             action='write', iostat=ierr)
        
        if (ierr /= 0) then
            if (this%verbose) then
                print *, "[ERROR] Cannot open file: ", trim(filename)
            end if
            return
        end if
        
        ! 写入头部信息（类似Julia的输出格式）
        write(unit, '(A)') "========================================"
        write(unit, '(A)') "CFD SOLVER RESULTS (Fortran)"
        write(unit, '(A)') "========================================"
        write(unit, '(A, A)') "Solver: ", trim(results%solver_name)
        write(unit, '(A, A)') "Scheme: ", trim(results%scheme)
        write(unit, '(A, I0)') "Order: ", results%order
        write(unit, '(A, I0)') "RK Order: ", results%rk_order
        write(unit, '(A, ES15.8)') "Final Time: ", results%final_time
        write(unit, '(A, ES15.8)') "Current Time: ", results%current_time
        write(unit, '(A, I0)') "Total Steps: ", results%total_steps
        write(unit, '(A, I0)') "Solver State: ", results%solver_state
        write(unit, '(A, I0)') "Grid Points: ", n
        write(unit, '(A)') "========================================"
        write(unit, '(A)') "DATA: x, numerical, analytical"
        write(unit, '(A)') "========================================"
        
        ! 写入数据
        do i = 1, n
            write(unit, '(3ES20.12)') results%x(i), results%numerical(i), results%analytical(i)
        end do
        
        ! 关闭文件
        close(unit)
        
        if (this%verbose) then
            print *, "[RESULTS] Saved ", n, " data points to ", trim(filename)
        end if
    end subroutine results_saver_save_text
    
    ! 二进制保存（可选）
    subroutine results_saver_save_binary(this, results, filename)
        class(results_saver), intent(in) :: this
        type(cfd_results), intent(in) :: results
        character(len=*), intent(in) :: filename
        
        ! 暂时实现文本格式，二进制格式可后续添加
        if (this%verbose) then
            print *, "[INFO] Binary save not implemented, using text format"
        end if
        call this%save_text(results, filename)
    end subroutine results_saver_save_binary
    
    ! 加载结果（暂时简单实现）
    subroutine results_saver_load(this, filename, results)
        class(results_saver), intent(in) :: this
        character(len=*), intent(in) :: filename
        type(cfd_results), intent(out) :: results
        
        ! 简化：只打印文件信息
        if (this%verbose) then
            print *, "[RESULTS] Would load from: ", trim(filename)
            print *, "  Note: Load functionality needs implementation"
        end if
        
        ! 初始化结果结构以避免未初始化警告
        results%solver_name = ""
        results%scheme = ""
        results%order = 0
        results%rk_order = 0
        results%final_time = 0.0_wp
        results%current_time = 0.0_wp
        results%total_steps = 0
        results%solver_state = 0
    end subroutine results_saver_load
    
end module results_module
