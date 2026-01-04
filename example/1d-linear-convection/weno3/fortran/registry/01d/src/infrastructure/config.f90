! src/infrastructure/config.f90
module config_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private
    public :: wp, cfd_config, config_print, config_with_reconstruction
    
    ! CFD configuration type
    type :: cfd_config
        character(len=20) :: ic_type = "step"
        character(len=20) :: recon_scheme = "eno"
        character(len=20) :: flux_type = "rusanov"
        integer :: rk_order = 1
        real(wp) :: wave_speed = 1.0_wp
        real(wp) :: final_time = 0.625_wp
        real(wp) :: dt = 0.025_wp
        character(len=20) :: boundary_type = "periodic"
        real(wp) :: left_boundary_value = 1.0_wp
        real(wp) :: right_boundary_value = 2.0_wp
        integer :: spatial_order = 2
        logical :: verbose = .true.
    end type cfd_config
    
    ! Interfaces
    interface config_print
        module procedure config_print_proc
    end interface
    
    interface config_with_reconstruction
        module procedure config_with_reconstruction_proc
    end interface
    
contains
    
    subroutine config_print_proc(this)
        type(cfd_config), intent(in) :: this
        
        print *, "=== CFD Configuration ==="
        print *, "Initial condition: ", trim(this%ic_type)
        print *, "Reconstruction: ", trim(this%recon_scheme), " (order:", this%spatial_order, ")"
        print *, "Flux type: ", trim(this%flux_type)
        print *, "Time integration: RK", this%rk_order
        print *, "Wave speed: ", this%wave_speed
        print *, "Final time: ", this%final_time
        print *, "Time step: ", this%dt
        print *, "Boundary: ", trim(this%boundary_type)
        if (trim(this%boundary_type) == 'dirichlet') then
            print *, "  Dirichlet values: [", this%left_boundary_value, ", ", &
                     this%right_boundary_value, "]"
        end if
        print *, "==============================="
    end subroutine config_print_proc
    
    subroutine config_with_reconstruction_proc(this, scheme, order)
        type(cfd_config), intent(inout) :: this
        character(len=*), intent(in) :: scheme
        integer, optional, intent(in) :: order
        
        character(len=20) :: scheme_lower
        
        ! Convert to lowercase
        scheme_lower = scheme
        call to_lower_inplace(scheme_lower)
        this%recon_scheme = trim(adjustl(scheme_lower))
        
        ! Set order
        if (present(order)) then
            this%spatial_order = order
        else
            ! Smart defaults
            if (index(this%recon_scheme, 'weno') > 0) then
                this%spatial_order = 5
            else if (trim(this%recon_scheme) == 'eno') then
                this%spatial_order = 3
            else
                error stop "[ERROR] Unsupported reconstruction scheme: " // trim(this%recon_scheme)
            end if
        end if
        
        if (this%verbose) then
            print *, "[CONFIG] Reconstruction: ", trim(this%recon_scheme), &
                     " Order: ", this%spatial_order
        end if
    end subroutine config_with_reconstruction_proc
    
    subroutine to_lower_inplace(str)
        character(len=*), intent(inout) :: str
        integer :: i
        
        do i = 1, len_trim(str)
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
                str(i:i) = char(ichar(str(i:i)) + 32)
            end if
        end do
    end subroutine to_lower_inplace
    
end module config_module