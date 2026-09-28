! src/boundary/factory.f90
module boundary_factory_module
    use base_modules, only: wp, ip
    use boundary_base_module, only: boundary_condition
    use periodic_boundary_module, only: periodic_boundary
    implicit none
    private
    
    public :: create_boundary_condition, initialize_boundary_factory
    
contains

    subroutine initialize_boundary_factory()
        print *, "[BOUNDARY FACTORY] Boundary conditions placeholder"
    end subroutine
    
    subroutine create_boundary_condition(bc_type, bc_instance, left_val, right_val)
        character(len=*), intent(in) :: bc_type
        class(boundary_condition), allocatable, intent(out) :: bc_instance
        real(wp), optional, intent(in) :: left_val, right_val
        
        ! 暂时只实现周期性边界
        allocate(periodic_boundary :: bc_instance)
        
        select case (trim(bc_type))
        case ("periodic")
            select type(bc => bc_instance)
            type is (periodic_boundary)
                bc%name = "periodic"
            end select
            
        case default
            print *, "[WARNING] Using periodic as default boundary"
            select type(bc => bc_instance)
            type is (periodic_boundary)
                bc%name = "periodic"
            end select
        end select
    end subroutine
    
end module boundary_factory_module