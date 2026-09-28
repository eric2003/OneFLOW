! src/boundary/boundary_base.f90 (修正版)
module boundary_base_module
    use base_modules, only: wp, ip
    implicit none
    private
    
    type, abstract, public :: boundary_condition
        character(len=:), allocatable :: name
    contains
        procedure(apply_interface), deferred :: apply
        procedure :: get_name => bc_get_name
    end type
    
    abstract interface
        subroutine apply_interface(this, u, nghosts, ist, ied)
            import :: boundary_condition, wp, ip
            class(boundary_condition), intent(in) :: this
            real(wp), intent(inout) :: u(:)
            integer(ip), intent(in) :: nghosts, ist, ied
        end subroutine apply_interface
    end interface
    
contains

    function bc_get_name(this) result(name)
        class(boundary_condition), intent(in) :: this
        character(len=:), allocatable :: name
        if (allocated(this%name)) then
            name = this%name
        else
            name = "unnamed"
        end if
    end function
    
end module boundary_base_module