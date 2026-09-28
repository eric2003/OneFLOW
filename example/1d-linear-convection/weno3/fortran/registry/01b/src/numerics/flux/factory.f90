! src/numerics/flux/factory.f90
module flux_factory_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use registry_module, only: create_component, has_component
    use flux_base_module, only: flux_calculator_base
    implicit none
    
    private
    public :: wp, flux_factory_create
    
contains
    
    subroutine flux_factory_create(flux_name, flux_calculator)
        character(len=*), intent(in) :: flux_name
        class(flux_calculator_base), allocatable, intent(out) :: flux_calculator
        
        class(*), allocatable :: instance
        character(len=20) :: name_lower
        integer :: i
        
        ! Convert to lowercase
        name_lower = flux_name
        do i = 1, len_trim(name_lower)
            if (name_lower(i:i) >= 'A' .and. name_lower(i:i) <= 'Z') then
                name_lower(i:i) = char(ichar(name_lower(i:i)) + 32)
            end if
        end do
        
        ! Check if registered
        if (.not. has_component("flux", trim(name_lower))) then
            print *, "[ERROR] Flux calculator not registered: ", trim(flux_name)
            return
        end if
        
        ! Create instance
        call create_component("flux", trim(name_lower), instance)
        
        ! Type conversion
        select type (inst => instance)
        type is (flux_calculator_base)
            allocate(flux_calculator, source=inst)
        class default
            error stop "[ERROR] Created flux calculator has wrong type"
        end select
    end subroutine flux_factory_create
    
end module flux_factory_module