! src/infrastructure/config.f90
module config_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use registry_module, only: register_component
    implicit none
    
    private
    public :: wp, cfd_config
    
    ! CFD配置类型
    type :: cfd_config
        character(len=20) :: ic_type = "step"
        character(len=20) :: recon_scheme = "eno"
        character(len=20) :: flux_type = "rusanov"
        integer :: rk_order = 1
        real(wp) :: wave_speed = 1.0_wp
        real(wp) :: final_time = 0.625_wp
        real(wp) :: dt = 0.025_wp
        character(len=20) :: boundary_type = "periodic"
        real(wp) :: left_boundary_value = 1.0_wp
        real(wp) :: right_boundary_value = 2.0_wp
        integer :: spatial_order = 2
        logical :: verbose = .true.
    contains
        procedure :: with_reconstruction => config_with_reconstruction
        procedure :: with_boundary => config_with_boundary
        procedure :: print => config_print
    end type cfd_config
    
contains
    
    subroutine config_with_reconstruction(this, scheme, order)
        class(cfd_config), intent(inout) :: this
        character(len=*), intent(in) :: scheme
        integer, optional, intent(in) :: order
        
        character(len=20) :: scheme_lower
        
        ! 转换为小写
        scheme_lower = scheme
        call to_lower_inplace(scheme_lower)
        this%recon_scheme = trim(adjustl(scheme_lower))
        
        ! 设置阶数
        if (present(order)) then
            this%spatial_order = order
        else
            ! 智能默认
            if (index(this%recon_scheme, 'weno') > 0) then
                this%spatial_order = 5
            else if (trim(this%recon_scheme) == 'eno') then
                this%spatial_order = 3
            else
                error stop "[ERROR] 不支持的重建格式: " // trim(this%recon_scheme)
            end if
        end if
        
        if (this%verbose) then
            print *, "[CONFIG] 重建方案: ", trim(this%recon_scheme), &
                     " 阶数: ", this%spatial_order
        end if
    end subroutine config_with_reconstruction
    
    subroutine config_with_boundary(this, bc_type, left_value, right_value)
        class(cfd_config), intent(inout) :: this
        character(len=*), intent(in) :: bc_type
        real(wp), optional, intent(in) :: left_value
        real(wp), optional, intent(in) :: right_value
        
        this%boundary_type = trim(adjustl(bc_type))
        
        if (present(left_value)) then
            this%left_boundary_value = left_value
        end if
        
        if (present(right_value)) then
            this%right_boundary_value = right_value
        end if
        
        if (this%verbose) then
            print *, "[CONFIG] 边界条件: ", trim(this%boundary_type), &
                     " 左值: ", this%left_boundary_value, &
                     " 右值: ", this%right_boundary_value
        end if
    end subroutine config_with_boundary
    
    subroutine config_print(this)
        class(cfd_config), intent(in) :: this
        
        print *, "=== CFD 配置 ==="
        print *, "初始条件: ", trim(this%ic_type)
        print *, "重建方案: ", trim(this%recon_scheme), " (order:", this%spatial_order, ")"
        print *, "通量类型: ", trim(this%flux_type)
        print *, "时间积分: RK", this%rk_order
        print *, "波速: ", this%wave_speed
        print *, "最终时间: ", this%final_time
        print *, "时间步长: ", this%dt
        print *, "边界条件: ", trim(this%boundary_type)
        if (trim(this%boundary_type) == 'dirichlet') then
            print *, "  Dirichlet值: [", this%left_boundary_value, ", ", &
                     this%right_boundary_value, "]"
        end if
        print *, "=============================="
    end subroutine config_print
    
    subroutine to_lower_inplace(str)
        character(len=*), intent(inout) :: str
        integer :: i
        
        do i = 1, len_trim(str)
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
                str(i:i) = char(ichar(str(i:i)) + 32)
            end if
        end do
    end subroutine to_lower_inplace
    
end module config_module