! src/core/factory_interfaces.f90
module factory_interfaces
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    private
    public :: wp, factory_procedure

    ! Factory procedure interface
    abstract interface
        subroutine factory_procedure(instance)
            import :: wp
            class(*), allocatable, intent(out) :: instance
        end subroutine factory_procedure
    end interface

end module factory_interfaces