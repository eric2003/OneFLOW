! ==================== residual.f90 ====================
! Residual computation module for OneFLOW-CFD

module residual_module
    use kinds, only: dp
    use config_module, only: CfdConfigType
    use domain_module, only: ComputationalDomainType
    use solution_module, only: SolutionType
    use flux_module, only: inviscid_flux
    use recon_module, only: ReconType, reconstruction
    implicit none
    
    private
    
    public :: residual

contains

    ! Compute residual (flux divergence)
    subroutine residual(recon, q, config, domain, solution)
		type(ReconType), intent(inout) :: recon
        real(dp), intent(in) :: q(:)
		type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        type(SolutionType), intent(inout) :: solution
        
        integer :: i
        
        ! Reconstruction using current solution
        call reconstruction(recon, q, config, domain, solution)
        
        ! Compute fluxes
        call inviscid_flux(solution%q_face_left, solution%q_face_right, &
                          solution%flux, config, domain%mesh)
        
        ! Compute residual
        do i = 1, domain%mesh%ncells
            solution%res(i) = -(solution%flux(i+1) - solution%flux(i)) / &
                                   domain%mesh%dx
        end do
    end subroutine residual
    
end module residual_module