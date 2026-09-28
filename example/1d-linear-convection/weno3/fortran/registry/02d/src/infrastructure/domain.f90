! src/infrastructure/domain.f90
module domain_module
    use base_modules, only: wp, ip
    use mesh_module, only: mesh_type
    use config_module, only: cfd_config
    
    implicit none
    private
    public :: wp, ip, domain_type, domain_create, is_physical_cell
    
    type :: domain_type
        type(cfd_config), pointer :: config => null()
        type(mesh_type), pointer :: mesh => null()
        integer(ip) :: nghosts = 0
        integer(ip) :: ist = 1      ! 1-based start index
        integer(ip) :: ied = 1      ! 1-based end index (exclusive)
        integer(ip) :: ntcells = 0
    contains
        procedure :: print_info => domain_print_info
    end type domain_type
    
contains

    function domain_create(config, mesh) result(domain)
        type(cfd_config), target, intent(in) :: config
        type(mesh_type), target, intent(in) :: mesh
        type(domain_type) :: domain
        
        domain%config => config
        domain%mesh => mesh
        
        ! 计算ghost层数（参考Julia的_calc_nghosts）
        domain%nghosts = calc_nghosts(config)
        domain%ist = domain%nghosts + 1
        domain%ied = domain%ist + mesh%ncells - 1
        domain%ntcells = mesh%ncells + 2 * domain%nghosts
        
        if (config%verbose) then
            print *, "[DOMAIN] Created with nghosts = ", domain%nghosts
            print *, "         Physical cells: ", domain%ist, " to ", domain%ied
            print *, "         Total cells: ", domain%ntcells
        end if
    end function domain_create
    
    function calc_nghosts(config) result(nghosts)
        type(cfd_config), intent(in) :: config
        integer(ip) :: nghosts
        
        character(len=max_name_length) :: scheme
        
        scheme = trim(config%recon_scheme)
        
        if (scheme == "eno") then
            nghosts = config%spatial_order
        else if (index(scheme, 'weno') > 0) then
            nghosts = config%spatial_order / 2 + 1
        else
            print *, "[ERROR] Unknown reconstruction scheme: ", trim(scheme)
            nghosts = 2  ! 默认值
        end if
    end function calc_nghosts
    
    logical function is_physical_cell(this, idx)
        class(domain_type), intent(in) :: this
        integer(ip), intent(in) :: idx
        is_physical_cell = (idx >= this%ist .and. idx < this%ied)
    end function is_physical_cell
    
    subroutine domain_print_info(this)
        class(domain_type), intent(in) :: this
        
        print *, "=== Domain Information ==="
        print *, "Ghost layers: ", this%nghosts
        print *, "Physical cells: ", this%ist, " to ", this%ied - 1
        print *, "Total cells: ", this%ntcells
        print *, "Mesh cells: ", this%mesh%ncells
        print *, "=========================="
    end subroutine domain_print_info

end module domain_module