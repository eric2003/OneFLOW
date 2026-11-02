! src/solver/residual.f90 (增强版)
module residual_module
    use base_modules, only: wp, ip
    use domain_module, only: domain_type
    use solution_module, only: solution_type
    use mesh_module, only: mesh_type
    use config_module, only: cfd_config
    use reconstructor_base_module, only: reconstructor_base
    use flux_base_module, only: flux_calculator_base
    
    implicit none
    private
    
    type, public :: residual_calculator
        ! 配置和域
        type(cfd_config), pointer :: config => null()
        type(mesh_type), pointer :: mesh => null()
        type(domain_type), pointer :: domain => null()
        
        ! 数值组件
        class(reconstructor_base), pointer :: reconstructor => null()
        class(flux_calculator_base), pointer :: flux_calc => null()
        
        ! 工作数组
        real(wp), pointer :: qL(:) => null()
        real(wp), pointer :: qR(:) => null()
        real(wp), pointer :: flux(:) => null()
        real(wp), pointer :: res(:) => null()
        
        real(wp) :: dx = 0.0_wp
        logical :: initialized = .false.
        
    contains
        procedure :: initialize => residual_init
        procedure :: compute => residual_compute
        procedure :: cleanup => residual_cleanup
    end type residual_calculator
    
contains

    subroutine residual_init(this, config, mesh, domain, reconstructor, flux_calc, &
                            solution)
        class(residual_calculator), intent(inout) :: this
        type(cfd_config), target, intent(in) :: config
        type(mesh_type), target, intent(in) :: mesh
        type(domain_type), target, intent(in) :: domain
        class(reconstructor_base), target, intent(in) :: reconstructor
        class(flux_calculator_base), target, intent(in) :: flux_calc
        type(solution_type), target, intent(in) :: solution
        
        this%config => config
        this%mesh => mesh
        this%domain => domain
        this%reconstructor => reconstructor
        this%flux_calc => flux_calc
        
        ! 链接工作数组
        this%qL => solution%q_face_left
        this%qR => solution%q_face_right
        this%flux => solution%flux
        this%res => solution%res
        
        ! 计算网格间距
        this%dx = mesh%dx
        
        this%initialized = .true.
        
        if (config%verbose) then
            print *, "[RESIDUAL] Initialized residual calculator"
            print *, "  Scheme: ", trim(reconstructor%name)
            print *, "  Flux: ", trim(flux_calc%name)
            print *, "  dx: ", this%dx
        end if
    end subroutine residual_init
    
    subroutine residual_compute(this, u)
        class(residual_calculator), intent(inout) :: this
        real(wp), intent(in) :: u(:)
        
        integer :: i, n_faces, n_cells
        real(wp) :: f_left, f_right
        
        if (.not. this%initialized) then
            if (this%config%verbose) then
                print *, "[ERROR] Residual calculator not initialized"
            end if
            return
        end if
        
        n_faces = this%mesh%nnodes
        n_cells = this%mesh%ncells
        
        ! 阶段1: 重构界面值（简化实现）
        ! 在实际实现中，这里应该调用重构器
        if (this%config%verbose .and. mod(this%domain%current_step, 100) == 0) then
            print *, "  [RESIDUAL] Reconstructing with ", trim(this%reconstructor%name)
        end if
        
        ! 简化重构：线性插值
        do i = 1, n_faces - 1
            this%qL(i+1) = 0.5_wp * (u(this%domain%ist + i - 2) + &
                                     u(this%domain%ist + i - 1))
            this%qR(i) = this%qL(i+1)
        end do
        
        ! 周期边界处理
        this%qL(1) = 0.5_wp * (u(this%domain%ied - 1) + u(this%domain%ist))
        this%qR(n_faces) = this%qL(1)
        
        ! 阶段2: 计算通量（简化实现）
        ! 在实际实现中，这里应该调用通量计算器
        if (this%config%verbose .and. mod(this%domain%current_step, 100) == 0) then
            print *, "  [RESIDUAL] Computing fluxes with ", trim(this%flux_calc%name)
        end if
        
        ! 简化通量：Lax-Friedrichs
        do i = 1, n_faces
            f_left = this%config%wave_speed * this%qL(i)
            f_right = this%config%wave_speed * this%qR(i)
            this%flux(i) = 0.5_wp * (f_left + f_right) - &
                          0.5_wp * abs(this%config%wave_speed) * (this%qR(i) - this%qL(i))
        end do
        
        ! 阶段3: 计算残差
        do i = 1, n_cells
            this%res(i) = -(this%flux(i+1) - this%flux(i)) / this%dx
        end do
        
        ! 调试输出
        if (this%config%verbose .and. mod(this%domain%current_step, 100) == 0) then
            print *, "  [RESIDUAL] Residual range: ", &
                    minval(this%res), " to ", maxval(this%res)
        end if
    end subroutine residual_compute
    
    subroutine residual_cleanup(this)
        class(residual_calculator), intent(inout) :: this
        this%initialized = .false.
        
        ! 清空指针
        this%config => null()
        this%mesh => null()
        this%domain => null()
        this%reconstructor => null()
        this%flux_calc => null()
        this%qL => null()
        this%qR => null()
        this%flux => null()
        this%res => null()
    end subroutine residual_cleanup

end module residual_module