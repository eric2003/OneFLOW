! src/core/registry_initializer.f90 (新增文件)
module registry_initializer_module
    use registry_module, only: register_component_simple
    implicit none
    private
    public :: initialize_default_registry
    
contains

    subroutine initialize_default_registry()
        ! 注册重构器
        call register_component_simple("reconstructor", "eno", order=3)
        call register_component_simple("reconstructor", "weno3", order=3)
        call register_component_simple("reconstructor", "weno5", order=5)
        
        ! 注册通量计算器
        call register_component_simple("flux", "rusanov")
        call register_component_simple("flux", "engquist-osher")
        
        ! 注册边界条件
        call register_component_simple("boundary", "periodic")
        call register_component_simple("boundary", "dirichlet")
        call register_component_simple("boundary", "neumann")
        
        ! 注册时间积分器
        call register_component_simple("integrator", "rk1", order=1)
        call register_component_simple("integrator", "rk2", order=2)
        call register_component_simple("integrator", "rk3", order=3)
        
        ! 注册方程
        call register_component_simple("equation", "linear_advection")
        
        ! 注册问题
        call register_component_simple("problem", "linear_advection")
        
        print *, "[REGISTRY] Default components registered"
    end subroutine initialize_default_registry

end module registry_initializer_module