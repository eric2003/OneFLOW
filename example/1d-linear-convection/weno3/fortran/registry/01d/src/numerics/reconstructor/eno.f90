! src/numerics/reconstructor/eno.f90
module eno_reconstructor_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use reconstructor_base_module, only: reconstructor_base, wp
    use registry_module, only: register_component_with_factory
    implicit none
    
    private
    public :: wp, eno_reconstructor, create_eno
    
    type, extends(reconstructor_base) :: eno_reconstructor
        real(wp) :: epsilon = 1.0e-6_wp
    contains
        procedure :: reconstruct => eno_reconstruct
        procedure :: info => eno_info
    end type eno_reconstructor
    
contains
    
    ! Factory function
    subroutine create_eno(instance)
        class(*), allocatable, intent(out) :: instance
        type(eno_reconstructor), allocatable :: eno
        
        allocate(eno)
        eno%name = "ENO"
        eno%order = 3
        eno%epsilon = 1.0e-6_wp
        
        call move_alloc(eno, instance)
    end subroutine create_eno
    
    ! ENO reconstruction (simplified 3rd order)
    subroutine eno_reconstruct(this, q, qL, qR)
        class(eno_reconstructor), intent(in) :: this
        real(wp), intent(in) :: q(:)
        real(wp), intent(out) :: qL(:)
        real(wp), intent(out) :: qR(:)
        
        integer :: i, n
        
        n = size(q) - 1  ! Number of interfaces
        
        if (n /= size(qL) .or. n /= size(qR)) then
            print *, "[ERROR] Array size mismatch in eno_reconstruct"
            return
        end if
        
        ! Simple implementation: 3rd order ENO
        do i = 1, n
            if (i >= 2 .and. i <= n-1) then
                ! 3-point stencil for 3rd order
                qL(i) = (1.0_wp/3.0_wp)*q(i-1) - (7.0_wp/6.0_wp)*q(i) + (11.0_wp/6.0_wp)*q(i+1)
                qR(i) = (-1.0_wp/6.0_wp)*q(i) + (5.0_wp/6.0_wp)*q(i+1) + (1.0_wp/3.0_wp)*q(i+2)
            else
                ! Use simple upwind near boundaries
                qL(i) = q(i)
                qR(i) = q(i+1)
            end if
        end do
    end subroutine eno_reconstruct
    
    subroutine eno_info(this)
        class(eno_reconstructor), intent(in) :: this
        call reconstructor_info(this)
        print *, "  Type: ENO reconstructor"
        print *, "  Epsilon: ", this%epsilon
    end subroutine eno_info
    
end module eno_reconstructor_module