! src/core/factory_integrated.f90
module factory_integrated
    use base_modules, only: wp, ip
    use registry_module, only: register_component_simple
    implicit none
    private
    
    public :: register_all_components
    
contains

    subroutine register_all_components()
        ! 方程
        call register_component_simple("equation", "linear_advection")
        
        ! 问题
        call register_component_simple("problem", "linear_advection")
        
        ! 重构器
        call register_component_simple("reconstructor", "eno")
        call register_component_simple("reconstructor", "weno3")
        call register_component_simple("reconstructor", "weno5")
        
        ! 通量
        call register_component_simple("flux", "rusanov")
        call register_component_simple("flux", "engquist_osher")
        
        ! 边界条件
        call register_component_simple("boundary", "periodic")
        call register_component_simple("boundary", "dirichlet")
        call register_component_simple("boundary", "neumann")
        
        ! 初始条件
        call register_component_simple("initial_condition", "step")
        call register_component_simple("initial_condition", "sin")
        call register_component_simple("initial_condition", "gaussian")
        
        print *, "[FACTORY] All components registered"
    end subroutine
    
end module