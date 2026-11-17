! src/numerics/time_integration/factory.f90
module time_integrator_factory_module
    use base_modules, only: wp, ip
    use time_integration_base_module, only: time_integrator_base
    use rk1_integrator_module, only: rk1_integrator
    use rk2_integrator_module, only: rk2_integrator
    use rk3_integrator_module, only: rk3_integrator
    
    implicit none
    private
    
    public :: create_time_integrator
    
contains

    function create_time_integrator(rk_order, config, domain, solution, residual_calc, bc) &
            result(integrator)
        integer, intent(in) :: rk_order
        type(cfd_config), target, intent(in) :: config
        type(domain_type), target, intent(in) :: domain
        type(solution_type), target, intent(in) :: solution
        type(residual_calculator), target, intent(in) :: residual_calc
        class(boundary_condition), target, intent(in) :: bc
        class(time_integrator_base), allocatable :: integrator
        
        select case (rk_order)
        case (1)
            allocate(rk1_integrator :: integrator)
            
        case (2)
            allocate(rk2_integrator :: integrator)
            
        case (3)
            allocate(rk3_integrator :: integrator)
            
        case default
            if (config%verbose) then
                print *, "[WARNING] Unsupported RK order: ", rk_order, ", using RK2 as default"
            end if
            allocate(rk2_integrator :: integrator)
        end select
        
        ! 初始化积分器
        call integrator%initialize(config, domain, solution, residual_calc, bc)
        
        if (config%verbose) then
            print *, "[TIME INTEGRATION] Created RK", rk_order, " integrator"
        end if
    end function create_time_integrator

end module time_integrator_factory_module