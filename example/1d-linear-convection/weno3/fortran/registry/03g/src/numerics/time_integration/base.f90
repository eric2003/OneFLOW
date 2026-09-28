! src/numerics/time_integration/base.f90
module time_integration_base_module
    use base_modules, only: wp, ip
    use config_module, only: cfd_config
    use domain_module, only: domain_type
    use solution_module, only: solution_type
    use residual_module, only: residual_calculator
    use boundary_base_module, only: boundary_condition
    
    implicit none
    private
    
    type, abstract, public :: time_integrator_base
        type(cfd_config), pointer :: config => null()
        type(domain_type), pointer :: domain => null()
        type(solution_type), pointer :: solution => null()
        type(residual_calculator), pointer :: residual_calc => null()
        class(boundary_condition), pointer :: bc => null()
    contains
        procedure :: initialize => integrator_initialize
        procedure(step_interface), deferred :: step
        procedure :: compute_residual => integrator_compute_residual
        procedure :: apply_boundary => integrator_apply_boundary
        procedure :: map_idx => integrator_map_idx
    end type time_integrator_base
    
    abstract interface
        subroutine step_interface(this, dt)
            import :: time_integrator_base, wp
            class(time_integrator_base), intent(inout) :: this
            real(wp), intent(in) :: dt
        end subroutine step_interface
    end interface
    
contains

    subroutine integrator_initialize(this, config, domain, solution, residual_calc, bc)
        class(time_integrator_base), intent(inout) :: this
        type(cfd_config), target, intent(in) :: config
        type(domain_type), target, intent(in) :: domain
        type(solution_type), target, intent(in) :: solution
        type(residual_calculator), target, intent(in) :: residual_calc
        class(boundary_condition), target, intent(in) :: bc
        
        this%config => config
        this%domain => domain
        this%solution => solution
        this%residual_calc => residual_calc
        this%bc => bc
    end subroutine integrator_initialize
    
    subroutine integrator_compute_residual(this)
        class(time_integrator_base), intent(inout) :: this
        
        if (associated(this%residual_calc) .and. associated(this%solution)) then
            call this%residual_calc%compute(this%solution%u)
        end if
    end subroutine integrator_compute_residual
    
    subroutine integrator_apply_boundary(this)
        class(time_integrator_base), intent(inout) :: this
        
        if (associated(this%bc) .and. associated(this%solution)) then
            call this%bc%apply(this%solution%u, &
                               this%domain%nghosts, &
                               this%domain%ist, &
                               this%domain%ied - 1)
        end if
    end subroutine integrator_apply_boundary
    
    integer function integrator_map_idx(this, i) result(idx)
        class(time_integrator_base), intent(in) :: this
        integer(ip), intent(in) :: i
        idx = i - this%domain%ist + 1
    end function integrator_map_idx

end module time_integration_base_module