! src/manager/component_factory.f90 (完整文件)
module component_factory_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use config_module, only: cfd_config
    use reconstructor_base_module, only: reconstructor_base
    use flux_base_module, only: flux_calculator_base
    
    use eno_reconstructor_module, only: eno_reconstructor, create_eno_reconstructor
    use weno3_reconstructor_module, only: weno3_reconstructor, create_weno3_reconstructor
    use weno5_reconstructor_module, only: weno5_reconstructor, create_weno5_reconstructor
    use rusanov_flux_module, only: rusanov_flux, create_rusanov_flux
    
    implicit none
    private
    public :: wp, create_reconstructor, create_flux_calculator
    
    ! 错误代码
    integer, parameter :: CM_SUCCESS = 0
    integer, parameter :: CM_ERROR_UNKNOWN_SCHEME = 1
    integer, parameter :: CM_ERROR_UNKNOWN_FLUX = 2
    integer, parameter :: CM_ERROR_INVALID_ORDER = 3
    
contains

    ! ==================== 重构器创建 ====================
    
    function create_reconstructor(config, status) result(recon)
        type(cfd_config), intent(in) :: config
        integer, optional, intent(out) :: status
        class(reconstructor_base), allocatable :: recon
        
        character(len=20) :: scheme
        integer :: order, error_code
        
        scheme = trim(adjustl(config%recon_scheme))
        order = config%spatial_order
        
        error_code = CM_SUCCESS
        
        if (config%verbose) then
            print *, "[COMPONENT FACTORY] Creating reconstructor: ", scheme, " order=", order
        end if
        
        ! 处理"weno"作为WENO5的别名（与Julia一致）
        if (scheme == "weno" .and. order == 5) then
            scheme = "weno5"
        end if
        
        select case(scheme)
        case('eno')
            allocate(eno_reconstructor :: recon)
            select type(recon)
            type is(eno_reconstructor)
                recon = create_eno_reconstructor()
                recon%order = order  ! 覆盖默认阶数
            end select
            
        case('weno3')
            allocate(weno3_reconstructor :: recon)
            select type(recon)
            type is(weno3_reconstructor)
                recon = create_weno3_reconstructor()
                recon%order = order  ! 覆盖默认阶数
            end select
            
        case('weno5')
            allocate(weno5_reconstructor :: recon)
            select type(recon)
            type is(weno5_reconstructor)
                recon = create_weno5_reconstructor()
                recon%order = order  ! 覆盖默认阶数
            end select
            
        case default
            error_code = CM_ERROR_UNKNOWN_SCHEME
            if (config%verbose) then
                print *, "[ERROR] Unknown reconstructor scheme: ", scheme
                print *, "        Available: eno, weno3, weno5"
            end if
        end select
        
        ! 检查阶数有效性
        if (error_code == CM_SUCCESS) then
            if (order < 1) then
                error_code = CM_ERROR_INVALID_ORDER
                if (config%verbose) then
                    print *, "[ERROR] Invalid spatial order: ", order
                end if
            end if
        end if
        
        ! 设置状态码
        if (present(status)) then
            status = error_code
        else if (error_code /= CM_SUCCESS) then
            error stop "Reconstructor creation failed"
        end if
    end function create_reconstructor

    ! ==================== 通量计算器创建 ====================
    
    function create_flux_calculator(config, status) result(flux)
        type(cfd_config), intent(in) :: config
        integer, optional, intent(out) :: status
        class(flux_calculator_base), allocatable :: flux
        
        character(len=20) :: flux_type
        integer :: error_code
        
        flux_type = trim(adjustl(config%flux_type))
        error_code = CM_SUCCESS
        
        if (config%verbose) then
            print *, "[COMPONENT FACTORY] Creating flux calculator: ", flux_type
        end if
        
        select case(flux_type)
        case('rusanov')
            allocate(rusanov_flux :: flux)
            select type(flux)
            type is(rusanov_flux)
                flux = create_rusanov_flux()
                flux%wave_speed_default = config%wave_speed
            end select
            
        case default
            error_code = CM_ERROR_UNKNOWN_FLUX
            if (config%verbose) then
                print *, "[ERROR] Unknown flux type: ", flux_type
                print *, "        Available: rusanov"
            end if
        end select
        
        ! 设置状态码
        if (present(status)) then
            status = error_code
        else if (error_code /= CM_SUCCESS) then
            error stop "Flux calculator creation failed"
        end if
    end function create_flux_calculator
    
end module component_factory_module