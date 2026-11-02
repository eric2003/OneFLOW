! ==================== boundary.f90 ====================
! Boundary conditions module for OneFLOW-CFD

module boundary_module
    use kinds, only: dp
    implicit none
    
    private
    
    public :: periodic_boundary, boundary
    
contains
    
    ! ===================================================================
    ! Boundary Conditions
    ! ===================================================================
    
    ! Periodic boundary conditions
    subroutine periodic_boundary(u, nghosts, ist, ied)
        real(dp), intent(inout) :: u(:)
        integer, intent(in) :: nghosts, ist, ied
        integer :: ig
        
        ! Left ghost cells = right interior cells
        do ig = 1, nghosts
            u(ist-ig) = u(ied -ig)
        end do
        
        ! Right ghost cells = left interior cells
        do ig = 1, nghosts
            u(ied-1+ig) = u(ist-1+ig)
        end do
    end subroutine
    
    subroutine boundary(u, nghosts, ist, ied)
        real(dp), intent(inout) :: u(:)
        integer, intent(in) :: nghosts, ist, ied
        call periodic_boundary(u, nghosts, ist, ied)
    end subroutine
    
end module boundary_module