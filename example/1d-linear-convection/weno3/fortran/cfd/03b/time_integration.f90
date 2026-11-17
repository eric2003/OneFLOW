! ==================== time_integration.f90 ====================
! Flux computation module for OneFLOW-CFD

module time_integration_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use config_module, only: CfdConfigType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType, update_oldfield
    use recon_module, only: ReconType, reconstruction
    use boundary_module, only: boundary
    use residual_module, only: residual
    implicit none
    private
    
    public :: runge_kutta
    
contains
   ! ===================================================================
    ! Time Integration
    ! ===================================================================
    ! 1st-order Runge-Kutta (Euler)
    subroutine runge_kutta_1(recon, config, domain, solution)
		type(ReconType), intent(inout) :: recon
		type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        integer :: i, j
        real(dp) :: dt
        
        dt = config%dt
        
        ! Compute residual
        call residual(recon, solution%u, config, domain, solution)
        
        ! Update solution
        do i = domain%ist, domain%ied - 1
            j = i - domain%ist + 1
            solution%u(i) = solution%u(i) + dt * solution%res(j)
        end do
        
        ! Apply boundary conditions again
        call boundary(solution%u, domain%nghosts, domain%ist, domain%ied)
        
        ! Save old solution
        call update_oldfield(solution%un, solution%u, domain%ntcells)
    end subroutine runge_kutta_1
    
    ! 2nd-order Runge-Kutta (Heun)
    subroutine runge_kutta_2(recon, config, domain, solution)
		type(ReconType), intent(inout) :: recon
		type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        integer :: i, j
        real(dp) :: dt
        
        dt = config%dt
        
        ! Stage 1
        call residual(recon, solution%u, config, domain, solution)
        
        do i = domain%ist, domain%ied - 1
            j = i - domain%ist + 1
            solution%u(i) =solution%u(i) + dt * solution%res(j)
        end do
        call boundary(solution%u, domain%nghosts, domain%ist, domain%ied)
        
        ! Stage 2
        call residual(recon, solution%u, config, domain, solution)
        do i = domain%ist, domain%ied - 1
            j = i - domain%ist + 1
            solution%u(i) = 0.5_dp * solution%un(i) + &
                            0.5_dp * solution%u(i) + &
                            0.5_dp * dt * solution%res(j)
        end do
        call boundary(solution%u, domain%nghosts, domain%ist, domain%ied)
        
        call update_oldfield(solution%un, solution%u, domain%ntcells)
    end subroutine runge_kutta_2
    
    ! Runge-Kutta selection
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
        case default
            call runge_kutta_1(recon, config, domain, solution)
        end select
    end subroutine runge_kutta
    
end module time_integration_module