! src/numerics/time_integration/base_simple.f90
module time_integration_base_simple_module
    use base_modules, only: wp, ip
    use config_module, only: cfd_config
    use domain_module, only: domain_type
    use residual_simple_module, only: residual_calculator_simple
    
    implicit none
    private
    
    type, abstract, public :: time_integrator_base_simple
        type(cfd_config), pointer :: config => null()
        type(domain_type), pointer :: domain => null()
        type(residual_calculator_simple), pointer :: residual_calc => null()
    contains
        procedure :: initialize => integrator_simple_initialize
        procedure(step_interface), deferred :: step
        procedure :: compute_residual => integrator_simple_compute_residual
    end type time_integrator_base_simple
    
    abstract interface
        subroutine step_interface(this, u, dt)
            import :: time_integrator_base_simple, wp
            class(time_integrator_base_simple), intent(inout) :: this
            real(wp), intent(inout) :: u(:)
            real(wp), intent(in) :: dt
        end subroutine step_interface
    end interface
    
contains

    subroutine integrator_simple_initialize(this, config, domain, residual_calc)
        class(time_integrator_base_simple), intent(inout) :: this
        type(cfd_config), target, intent(in) :: config
        type(domain_type), target, intent(in) :: domain
        type(residual_calculator_simple), target, intent(in) :: residual_calc
        
        this%config => config
        this%domain => domain
        this%residual_calc => residual_calc
    end subroutine integrator_simple_initialize
    
    subroutine integrator_simple_compute_residual(this, u, res)
        class(time_integrator_base_simple), intent(inout) :: this
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: res(:)
        
        if (associated(this%residual_calc)) then
            call this%residual_calc%compute(u, res)
        end if
    end subroutine integrator_simple_compute_residual

end module time_integration_base_simple_module