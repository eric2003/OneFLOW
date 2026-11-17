! src/initial_condition/ic_base.f90
module ic_base_module
    use base_modules, only: wp, ip
    use solution_module, only: solution_type
    implicit none
    private
    
    type, abstract, public :: initial_condition
        character(len=:), allocatable :: name
    contains
        procedure :: get_name => ic_get_name
        procedure(evaluate_interface), deferred :: evaluate_at
        procedure(apply_interface), deferred :: apply
    end type
    
    abstract interface
        function evaluate_interface(this, x) result(u)
            import :: initial_condition, wp
            class(initial_condition), intent(in) :: this
            real(wp), intent(in) :: x(:)
            real(wp) :: u(size(x))
        end function evaluate_interface
        
        subroutine apply_interface(this, solution)
            import :: initial_condition, solution_type
            class(initial_condition), intent(in) :: this
            type(solution_type), intent(inout) :: solution
        end subroutine apply_interface
    end interface
    
contains

    function ic_get_name(this) result(name)
        class(initial_condition), intent(in) :: this
        character(len=:), allocatable :: name
        if (allocated(this%name)) then
            name = this%name
        else
            name = "unnamed"
        end if
    end function
    
end module ic_base_module