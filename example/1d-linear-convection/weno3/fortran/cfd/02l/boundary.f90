! ==================== boundary.f90 ====================
! Boundary conditions module for OneFLOW-CFD

module boundary_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
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
        integer :: i, j
        
        ! Left ghost cells = right interior cells
        do i = 1, nghosts
            j = ist - i
            u(j) = u(ied - nghosts + i - 1)
        end do
        
        ! Right ghost cells = left interior cells
        do i = 1, nghosts
            j = ied + i
            u(j) = u(ist + i - 1)
        end do
    end subroutine
    
    subroutine boundary(u, nghosts, ist, ied)
        real(dp), intent(inout) :: u(:)
        integer, intent(in) :: nghosts, ist, ied
        call periodic_boundary(u, nghosts, ist, ied)
    end subroutine
    
end module boundary_module