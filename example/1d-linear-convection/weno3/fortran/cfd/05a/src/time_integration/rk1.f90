! ==================== rk1.f90 ====================
! rk1 module for OneFLOW-CFD

module rk1_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType, update_oldfield
    use recon_module, only: ReconType
    use boundary_module, only: boundary
    use residual_module, only: residual	
    implicit none
    private
    public :: runge_kutta_1
    
contains
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
end module rk1_module