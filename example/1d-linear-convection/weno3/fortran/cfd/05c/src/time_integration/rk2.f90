! ==================== rk2.f90 ====================
! rk2 module for OneFLOW-CFD

module rk2_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType, update_oldfield
    use recon_module, only: ReconType
    use boundary_module, only: boundary
    use residual_module, only: residual
    implicit none
    private
    public :: runge_kutta_2
    
contains
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
end module rk2_module