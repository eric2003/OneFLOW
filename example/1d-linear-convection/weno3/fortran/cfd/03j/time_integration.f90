! ==================== time_integration.f90 ====================
! Time Integration computation module for OneFLOW-CFD

module time_integration_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    use recon_module, only: ReconType
    use boundary_module, only: boundary
    use residual_module, only: residual
    
    ! 导入各阶RK模块
    use rk1_module, only: runge_kutta_1
    use rk2_module, only: runge_kutta_2
    use rk3_module, only: runge_kutta_3	
	
    implicit none
    private
    public :: runge_kutta
    
contains
    subroutine runge_kutta(recon, config, domain, solution)
        type(ReconType), intent(inout) :: recon
        type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        select case(config%rk_order)
        case(1)
            call runge_kutta_1(recon, config, domain, solution)
        case(2)
            call runge_kutta_2(recon, config, domain, solution)
        case(3)
            call runge_kutta_3(recon, config, domain, solution)
        case default
            call runge_kutta_1(recon, config, domain, solution)
        end select
    end subroutine runge_kutta
end module time_integration_module    
