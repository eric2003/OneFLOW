! src/boundary/periodic.f90 (修正版)
module periodic_boundary_module
    use base_modules, only: wp, ip
    use boundary_base_module, only: boundary_condition
    implicit none
    private
    
    type, extends(boundary_condition), public :: periodic_boundary
    contains
        procedure :: apply => periodic_apply
        procedure :: get_name => bc_get_name
    end type
    
contains

    subroutine periodic_apply(this, u, nghosts, ist, ied)
        class(periodic_boundary), intent(in) :: this
        real(wp), intent(inout) :: u(:)
        integer(ip), intent(in) :: nghosts, ist, ied
        
        integer :: i
        
        ! 左ghost层：u[ist-1-i] = u[ied-1-i]
        do i = 0, nghosts-1
            u(ist-1-i) = u(ied-1-i)
        end do
        
        ! 右ghost层：u[ied+i] = u[ist+i]
        do i = 0, nghosts-1
            u(ied+i) = u(ist+i)
        end do
    end subroutine
    
    function bc_get_name(this) result(name)
        class(periodic_boundary), intent(in) :: this
        character(len=:), allocatable :: name
        if (allocated(this%name)) then
            name = this%name
        else
            name = "periodic"
        end if
    end function
    
end module periodic_boundary_module