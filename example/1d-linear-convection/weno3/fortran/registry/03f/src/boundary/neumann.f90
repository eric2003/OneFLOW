! src/boundary/neumann.f90
module neumann_boundary_module
    use base_modules, only: wp, ip
    use boundary_base_module, only: boundary_condition
    implicit none
    private
    
    type, extends(boundary_condition), public :: neumann_boundary
    contains
        procedure :: apply => neumann_apply
    end type
    
    interface neumann_boundary
        module procedure create_neumann_boundary
    end interface
    
contains

    type(neumann_boundary) function create_neumann_boundary() result(this)
        this%name = "neumann"
    end function
    
    subroutine neumann_apply(this, u, nghosts, ist, ied)
        class(neumann_boundary), intent(in) :: this
        real(wp), intent(inout) :: u(:)
        integer(ip), intent(in) :: nghosts, ist, ied
        
        integer :: i
        
        ! 左边界零梯度：u[ist-1-i] = u[ist+i]
        do i = 0, nghosts-1
            u(ist-1-i) = u(ist+i)
        end do
        
        ! 右边界零梯度：u[ied+i] = u[ied-1-i]
        do i = 0, nghosts-1
            u(ied+i) = u(ied-1-i)
        end do
    end subroutine
    
end module neumann_boundary_module