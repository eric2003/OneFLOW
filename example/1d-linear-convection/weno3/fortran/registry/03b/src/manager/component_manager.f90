! src/manager/component_manager.f90 (完整文件)
module component_manager_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use config_module, only: cfd_config
    use component_factory_module, only: create_reconstructor, create_flux_calculator
    use reconstructor_base_module, only: reconstructor_base
    use flux_base_module, only: flux_calculator_base
    
    implicit none
    private
    public :: wp, component_manager_info, validate_config
    public :: create_reconstructor, create_flux_calculator
    
contains

    ! ==================== 配置验证 ====================
    
    function validate_config(config) result(is_valid)
        type(cfd_config), intent(in) :: config
        logical :: is_valid
        
        integer :: status
        class(reconstructor_base), allocatable :: test_recon
        class(flux_calculator_base), allocatable :: test_flux
        
        is_valid = .false.
        
        ! 测试创建重构器
        test_recon = create_reconstructor(config, status)
        if (status /= 0) then
            if (config%verbose) then
                print *, "[CONFIG VALIDATION] Invalid reconstructor configuration"
            end if
            return
        end if
        
        ! 测试创建通量计算器
        test_flux = create_flux_calculator(config, status)
        if (status /= 0) then
            if (config%verbose) then
                print *, "[CONFIG VALIDATION] Invalid flux configuration"
            end if
            return
        end if
        
        ! 清理测试组件
        if (allocated(test_recon)) deallocate(test_recon)
        if (allocated(test_flux)) deallocate(test_flux)
        
        is_valid = .true.
        
        if (config%verbose) then
            print *, "[CONFIG VALIDATION] Configuration is valid"
        end if
    end function validate_config
    
    ! ==================== 信息显示 ====================
    
    subroutine component_manager_info()
        print *, "=== Component Manager ==="
        print *, "Available reconstructors:"
        print *, "  - eno   (orders: 1-7)"
        print *, "  - weno3 (order: 3)"
        print *, "  - weno5 (order: 5)"
        print *, ""
        print *, "Available flux calculators:"
        print *, "  - rusanov"
        print *, ""
        print *, "Features:"
        print *, "  - Configuration validation"
        print *, "  - Component creation from config"
        print *, "  - Error handling with status codes"
        print *, "========================="
    end subroutine component_manager_info
    
end module component_manager_module