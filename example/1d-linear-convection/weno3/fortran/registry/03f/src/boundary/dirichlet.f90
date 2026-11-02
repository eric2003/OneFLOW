! src/boundary/dirichlet.f90
module dirichlet_boundary_module
    use base_modules, only: wp, ip
    use boundary_base_module, only: boundary_condition
    implicit none
    private
    
    type, extends(boundary_condition), public :: dirichlet_boundary
        real(wp) :: left_value = 1.0_wp
        real(wp) :: right_value = 2.0_wp
    contains
        procedure :: apply => dirichlet_apply
        procedure :: set_values => dirichlet_set_values
    end type
    
    interface dirichlet_boundary
        module procedure create_dirichlet_boundary
    end interface
    
contains

    type(dirichlet_boundary) function create_dirichlet_boundary(left_val, right_val) result(this)
        real(wp), optional, intent(in) :: left_val, right_val
        
        this%name = "dirichlet"
        if (present(left_val)) this%left_value = left_val
        if (present(right_val)) this%right_value = right_val
    end function
    
    subroutine dirichlet_set_values(this, left_val, right_val)
        class(dirichlet_boundary), intent(inout) :: this
        real(wp), intent(in) :: left_val, right_val
        this%left_value = left_val
        this%right_value = right_val
    end subroutine
    
    subroutine dirichlet_apply(this, u, nghosts, ist, ied)
        class(dirichlet_boundary), intent(in) :: this
        real(wp), intent(inout) :: u(:)
        integer(ip), intent(in) :: nghosts, ist, ied
        
        integer :: i
        
        ! 左边界
        do i = 0, nghosts-1
            u(ist-1-i) = this%left_value
        end do
        
        ! 右边界
        do i = 0, nghosts-1
            u(ied+i) = this%right_value
        end do
    end subroutine
    
end module dirichlet_boundary_module