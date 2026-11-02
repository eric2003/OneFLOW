! OneFLOW-CFD Solver for 1D Convection Equation
! ENO/WENO Reconstruction Comparison - Single File Implementation

module cfd_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    private
    
    ! ===================================================================
    ! Forward Type Declarations
    ! ===================================================================
    type, public :: CfdConfigType
        character(len=10) :: recon_scheme = "eno"
        integer :: flux_type = 0
        integer :: rk_order = 1
        integer :: spatial_order = 3
        real(dp) :: wave_speed = 1.0_dp
        real(dp) :: final_time = 0.625_dp
        real(dp) :: dt = 0.025_dp  ! 减小时间步长
        real(dp) :: cfl = 0.5_dp   ! CFL数
    end type CfdConfigType
    
    type, public :: MeshType
        real(dp) :: xmin = 0.0_dp
        real(dp) :: xmax = 2.0_dp
        integer :: ncells = 40
        integer :: nnodes = 0
        integer :: nx = 0
        real(dp) :: L = 0.0_dp
        real(dp) :: dx = 0.0_dp
        real(dp), allocatable :: x(:), xcc(:)
    contains
        procedure :: init => mesh_init
    end type MeshType
    
    type, public :: ComputationalDomainType
        type(MeshType) :: mesh
        type(CfdConfigType) :: config
        integer :: nghosts = 0
        integer :: ist = 0
        integer :: ied = 0
        integer :: ntcells = 0
    contains
        procedure :: init => domain_init
    end type ComputationalDomainType
    
    type, public :: SolutionType
        type(ComputationalDomainType) :: domain
        real(dp), allocatable :: q_face_left(:), q_face_right(:)
        real(dp), allocatable :: flux(:), res(:)
        real(dp), allocatable :: u(:), un(:)
    contains
        procedure :: init => solution_init
    end type SolutionType
    
    ! Main CFD solver class
    type, public :: CfdType
        type(CfdConfigType) :: config
        type(ComputationalDomainType) :: domain
        type(SolutionType) :: solution
        class(*), allocatable :: reconstructor
    contains
        procedure :: init => cfd_init
    end type CfdType
    
    ! Abstract reconstructor base class
    type, abstract, public :: ReconstructorType
    contains
        procedure(reconstruct_interface), deferred, pass :: reconstruct
    end type ReconstructorType
    
    ! Define abstract interface
    abstract interface
        subroutine reconstruct_interface(this, q, cfd)
            import :: ReconstructorType, CfdType, dp
            class(ReconstructorType), intent(inout) :: this
            real(dp), intent(in) :: q(:)
            type(CfdType), intent(inout) :: cfd
        end subroutine reconstruct_interface
    end interface
    
    ! ENO reconstructor
    type, extends(ReconstructorType), public :: EnoReconstructorType
        integer :: spatial_order
        integer :: ntcells
        integer, allocatable :: lmc(:)
        real(dp), allocatable :: coef(:,:)
        real(dp), allocatable :: dd(:,:)
    contains
        procedure :: reconstruct => eno_reconstruct
        procedure :: init => eno_init
    end type EnoReconstructorType
    
    ! WENO reconstructor
    type, extends(ReconstructorType), public :: WenoReconstructorType
    contains
        procedure :: reconstruct => weno_reconstruct
    end type WenoReconstructorType
    
    ! ===================================================================
    ! Public Procedures
    ! ===================================================================
    public :: run_simulation, init_field, analytical_solution, performEnoWenoAnalysis
    
    ! ===================================================================
    ! Module Variables
    ! ===================================================================
    real(dp), parameter :: eps_weno = 1.0e-6_dp
    
contains
    
    ! ===================================================================
    ! Initialization Methods
    ! ===================================================================
    
    ! Mesh initialization
    subroutine mesh_init(this)
        class(MeshType), intent(inout) :: this
        integer :: i
        
        this%nnodes = this%ncells + 1
        this%nx = this%ncells
        this%L = this%xmax - this%xmin
        this%dx = this%L / real(this%ncells, dp)
        
        allocate(this%x(this%nnodes), this%xcc(this%ncells))
        
        ! Node coordinates
        do i = 1, this%nnodes
            this%x(i) = this%xmin + (i-1) * this%dx
        end do
        
        ! Cell center coordinates
        do i = 1, this%ncells
            this%xcc(i) = 0.5_dp * (this%x(i) + this%x(i+1))
        end do
    end subroutine mesh_init
    
    ! Domain initialization
    subroutine domain_init(this, mesh, config)
        class(ComputationalDomainType), intent(inout) :: this
        type(MeshType), intent(in) :: mesh
        type(CfdConfigType), intent(in) :: config
        
        this%mesh = mesh
        this%config = config
        
        ! Calculate ghost cells
        if (trim(config%recon_scheme) == "eno") then
            this%nghosts = config%spatial_order
        else if (trim(config%recon_scheme) == "weno") then
            this%nghosts = config%spatial_order / 2 + 1
        else
            error stop "Unknown reconstruction scheme"
        end if
        
        this%ist = this%nghosts + 1
        this%ied = this%ist + mesh%ncells - 1
        this%ntcells = mesh%ncells + 2 * this%nghosts
        
        print *, "Domain initialized:"
        print *, "  mesh.ncells = ", mesh%ncells
        print *, "  spatial_order = ", config%spatial_order
        print *, "  nghosts = ", this%nghosts
        print *, "  ist = ", this%ist, ", ied = ", this%ied
        print *, "  dx = ", mesh%dx
    end subroutine domain_init
    
    ! Solution initialization
    subroutine solution_init(this, domain)
        class(SolutionType), intent(inout) :: this
        type(ComputationalDomainType), intent(in) :: domain
        
        this%domain = domain
        
        allocate(this%q_face_left(domain%mesh%nnodes))
        allocate(this%q_face_right(domain%mesh%nnodes))
        allocate(this%flux(domain%mesh%nnodes))
        allocate(this%res(domain%mesh%ncells))
        allocate(this%u(domain%ntcells))
        allocate(this%un(domain%ntcells))
        
        this%q_face_left = 0.0_dp
        this%q_face_right = 0.0_dp
        this%flux = 0.0_dp
        this%res = 0.0_dp
        this%u = 0.0_dp
        this%un = 0.0_dp
    end subroutine solution_init
    
    ! ENO reconstructor initialization
    subroutine eno_init(this, spatial_order, ntcells)
        class(EnoReconstructorType), intent(inout) :: this
        integer, intent(in) :: spatial_order
        integer, intent(in) :: ntcells
        
        this%spatial_order = spatial_order
        this%ntcells = ntcells
        
        allocate(this%lmc(ntcells))
        allocate(this%coef(spatial_order+1, spatial_order))
        allocate(this%dd(spatial_order, ntcells))
        
        this%lmc = 0
        this%coef = 0.0_dp
        this%dd = 0.0_dp
        
        ! Initialize coefficients
        call init_coef(spatial_order, this%coef)
    end subroutine eno_init
    
    ! CFD solver initialization
    subroutine cfd_init(this, config, domain)
        class(CfdType), intent(inout) :: this
        type(CfdConfigType), intent(in) :: config
        type(ComputationalDomainType), intent(in) :: domain
        
        this%config = config
        this%domain = domain
        call this%solution%init(domain)
        
        ! Create reconstructor based on scheme
        if (trim(config%recon_scheme) == "eno") then
            allocate(EnoReconstructorType :: this%reconstructor)
            select type(rec => this%reconstructor)
            type is (EnoReconstructorType)
                call rec%init(config%spatial_order, domain%ntcells)
            end select
        else if (trim(config%recon_scheme) == "weno") then
            allocate(WenoReconstructorType :: this%reconstructor)
        else
            error stop "Unknown reconstruction scheme"
        end if
        
        ! Adjust time step based on CFL condition
        call calculate_dt(this)
    end subroutine cfd_init
    
    ! Calculate time step based on CFL condition
    subroutine calculate_dt(cfd)
        type(CfdType), intent(inout) :: cfd
        
        real(dp) :: dt_cfl
        
        ! CFL condition: dt <= CFL * dx / |wave_speed|
        dt_cfl = cfd%config%cfl * cfd%domain%mesh%dx / abs(cfd%config%wave_speed)
        
        if (cfd%config%dt > dt_cfl) then
            print *, "Adjusting time step for stability:"
            print *, "  Original dt = ", cfd%config%dt
            print *, "  CFL dt = ", dt_cfl
            cfd%config%dt = dt_cfl
            print *, "  Using dt = ", cfd%config%dt
        end if

    end subroutine calculate_dt
    
    ! ===================================================================
    ! Initial Conditions and Analytical Solution
    ! ===================================================================
    
    ! Initial condition: step function
    function initial_condition(x) result(u0)
        real(dp), intent(in) :: x
        real(dp) :: u0
        
        if (0.5_dp <= x .and. x <= 1.0_dp) then
            u0 = 2.0_dp
        else
            u0 = 1.0_dp
        end if
    end function initial_condition
    
    ! Analytical solution with periodic BC
    function analytical_solution(x, t, a, L) result(u)
        real(dp), intent(in) :: x, t, a, L
        real(dp) :: u, x_shifted
        
        x_shifted = mod(x - a * t + L, L)
        u = initial_condition(x_shifted)
    end function analytical_solution
    
    ! Initialize field with step function
    subroutine init_field(cfd)
        type(CfdType), intent(inout) :: cfd
        integer :: i, j
        
        do i = 1, cfd%domain%mesh%ncells
            if (0.5_dp <= cfd%domain%mesh%xcc(i) .and. cfd%domain%mesh%xcc(i) <= 1.0_dp) then
                cfd%solution%u(cfd%domain%ist + i - 1) = 2.0_dp
            else
                cfd%solution%u(cfd%domain%ist + i - 1) = 1.0_dp
            end if
        end do
        
        call boundary(cfd%solution%u, cfd)
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine init_field
    
    ! ===================================================================
    ! Boundary Conditions
    ! ===================================================================
    
    ! Periodic boundary conditions - CORRECTED VERSION
    subroutine periodic_boundary(u, cfd)
        real(dp), intent(inout) :: u(:)
        type(CfdType), intent(in) :: cfd
        integer :: i
        
        ! Copy interior cells to ghost cells
        ! Left ghost cells = right interior cells
        do i = 1, cfd%domain%nghosts
            u(cfd%domain%ist - i) = u(cfd%domain%ied - cfd%domain%nghosts + i - 1)
        end do
        
        ! Right ghost cells = left interior cells
        do i = 1, cfd%domain%nghosts
            u(cfd%domain%ied + i) = u(cfd%domain%ist + i - 1)
        end do
    end subroutine periodic_boundary
    
    ! Boundary condition wrapper
    subroutine boundary(u, cfd)
        real(dp), intent(inout) :: u(:)
        type(CfdType), intent(in) :: cfd
        
        call periodic_boundary(u, cfd)
    end subroutine boundary
    
    ! ===================================================================
    ! Reconstruction Methods
    ! ===================================================================
    
    ! Initialize ENO/WENO coefficients
    subroutine init_coef(spatial_order, coef)
        integer, intent(in) :: spatial_order
        real(dp), intent(out) :: coef(:,:)
        
        coef = 0.0_dp
        
        select case(spatial_order)
        case(1)
            coef(1,1) = 1.0_dp
            coef(2,1) = 1.0_dp
            
        case(2)
            coef(1,1:2) = [ 3.0_dp/2.0_dp, -1.0_dp/2.0_dp ]
            coef(2,1:2) = [ 1.0_dp/2.0_dp,  1.0_dp/2.0_dp ]
            coef(3,1:2) = [ -1.0_dp/2.0_dp, 3.0_dp/2.0_dp ]
            
        case(3)
            coef(1,1:3) = [ 11.0_dp/6.0_dp, -7.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
            coef(2,1:3) = [  1.0_dp/3.0_dp,  5.0_dp/6.0_dp, -1.0_dp/6.0_dp ]
            coef(3,1:3) = [ -1.0_dp/6.0_dp,  5.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
            coef(4,1:3) = [  1.0_dp/3.0_dp, -7.0_dp/6.0_dp, 11.0_dp/6.0_dp ]
            
        case default
            error stop "Unsupported spatial order"
        end select
    end subroutine init_coef
    
    subroutine first_order_reconstruct(this, q, cfd)
        class(EnoReconstructorType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i, j
        integer :: nghosts, ist, ied
        
        nghosts = cfd%domain%nghosts
        ist = cfd%domain%ist
        ied = cfd%domain%ied
        
        ! For now, use simple 2nd order reconstruction
        do i = ist, ied
            j = i - ist + 1  ! 1-based index for interfaces
            
            cfd%solution%q_face_left(j) = q(i-1)
            cfd%solution%q_face_right(j) = q(i)
        end do
    end subroutine first_order_reconstruct

    subroutine eno_reconstruct(this, q, cfd)
        class(EnoReconstructorType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(CfdType), intent(inout) :: cfd
		type(ComputationalDomainType), pointer :: domain
        
        integer :: i, j, m, k1, k2, r1, r2
        integer :: nghosts, ist, ied
        
        nghosts = cfd%domain%nghosts
        ist = cfd%domain%ist
        ied = cfd%domain%ied

		! 1. 差商计算 (dd[1,:] = q)
		this%dd(1, :) = q(:)
		do m = 2, this%spatial_order
			do j = 1, cfd%domain%ntcells - m + 1
				this%dd(m, j) = this%dd(m-1, j+1) - this%dd(m-1, j)
			end do
		end do
		
		! 2. 选择 smoothest stencil
		do i = ist - 1, ied   ! Python: range(ist-1, ied+1) → ied+1-1 = ied
			this%lmc(i) = i
			do m = 2, this%spatial_order
				if ( abs(this%dd(m, this%lmc(i) - 1) ) < abs(this%dd(m, this%lmc(i)))) then
					this%lmc(i) = this%lmc(i) - 1
                end if
            end do
        end do
		
		associate (  &
			q_face_left => cfd%solution%q_face_left,  &
			q_face_right => cfd%solution%q_face_right  &
		)  
			! 这里可以直接使用 q_face_left 和 q_face_right  
		! 3. 重构界面值
		do i = ist, ied
			j = i - ist + 1  ! 1-based index for interfaces
			k1 = this%lmc(i - 1)
			k2 = this%lmc(i)
			r1 = (i - 1) - k1 + 1
			r2 = i - k2 + 1
			q_face_left(j) = 0.0
			q_face_right(j) = 0.0
			do m = 1, this%spatial_order
				q_face_left(j) = q_face_left(j) + q(k1 + m - 1) * this%coef(r1 + 1, m)
				q_face_right(j) = q_face_right(j) + q(k2 + m - 1) * this%coef(r2, m)
			end do
		end do
		end associate
       
    end subroutine eno_reconstruct	
    
    ! WENO-3 nonlinear weights for left-biased stencil
    function wc3L(v1, v2, v3) result(f)
        real(dp), intent(in) :: v1, v2, v3
        real(dp) :: f, s0, s1, d0, d1, c0, c1, w0, w1, q0, q1
        
        ! Smoothness indicators
        s0 = (v3 - v2)**2
        s1 = (v2 - v1)**2
        
        ! Nonlinear weights
        d0 = 2.0_dp/3.0_dp
        d1 = 1.0_dp/3.0_dp
        
        c0 = d0 / ((eps_weno + s0)**2)
        c1 = d1 / ((eps_weno + s1)**2)
        
        w0 = c0 / (c0 + c1)
        w1 = c1 / (c0 + c1)
        
        ! Candidate stencils
        q0 = 0.5_dp * v2 + 0.5_dp * v3
        q1 = -0.5_dp * v1 + 1.5_dp * v2
        
        ! Reconstructed value
        f = w0 * q0 + w1 * q1
    end function wc3L
    
    ! WENO-3 nonlinear weights for right-biased stencil
    function wc3R(v1, v2, v3) result(f)
        real(dp), intent(in) :: v1, v2, v3
        real(dp) :: f, s0, s1, d0, d1, c0, c1, w0, w1, q0, q1
        
        ! Smoothness indicators
        s0 = (v2 - v1)**2
        s1 = (v3 - v2)**2
        
        ! Nonlinear weights
        d0 = 2.0_dp/3.0_dp
        d1 = 1.0_dp/3.0_dp
        
        c0 = d0 / ((eps_weno + s0)**2)
        c1 = d1 / ((eps_weno + s1)**2)
        
        w0 = c0 / (c0 + c1)
        w1 = c1 / (c0 + c1)
        
        ! Candidate stencils
        q0 = 0.5_dp * v1 + 0.5_dp * v2
        q1 = 1.5_dp * v2 - 0.5_dp * v3
        
        ! Reconstructed value
        f = w0 * q0 + w1 * q1
    end function wc3R
    
    ! WENO-3 reconstruction for left interface
    subroutine weno3L_periodic(cfd, u, qL)
        type(CfdType), intent(in) :: cfd
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: qL(:)
        
        integer :: i, j, nghosts, ist, ied
        
        nghosts = cfd%domain%nghosts
        ist = cfd%domain%ist
        ied = cfd%domain%ied
        
        do i = ist-1, ied-1
            j = i - (ist - 1) + 1
            qL(j) = wc3L(u(i-1), u(i), u(i+1))
        end do
	
    end subroutine weno3L_periodic
    
    ! WENO-3 reconstruction for right interface
    subroutine weno3R_periodic(cfd, u, qR)
        type(CfdType), intent(in) :: cfd
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: qR(:)
        
        integer :: i, j, nghosts, ist, ied
        
        nghosts = cfd%domain%nghosts
        ist = cfd%domain%ist
        ied = cfd%domain%ied
        
        do i = ist, ied
            j = i - ist + 1
            qR(j) = wc3R(u(i-1), u(i), u(i+1))
        end do

    end subroutine weno3R_periodic
    
    ! WENO reconstruction
    subroutine weno_reconstruct(this, q, cfd)
        class(WenoReconstructorType), intent(inout) :: this
        real(dp), intent(in) :: q(:)
        type(CfdType), intent(inout) :: cfd
        
        call weno3L_periodic(cfd, q, cfd%solution%q_face_left)
        call weno3R_periodic(cfd, q, cfd%solution%q_face_right)
    end subroutine weno_reconstruct
    
    ! General reconstruction wrapper
    subroutine reconstruction(q, cfd)
        real(dp), intent(in) :: q(:)
        type(CfdType), intent(inout) :: cfd
        
        select type(rec => cfd%reconstructor)
        type is (EnoReconstructorType)
            call rec%reconstruct(q, cfd)
        type is (WenoReconstructorType)
            call rec%reconstruct(q, cfd)
        class default
            error stop "Unknown reconstructor type"
        end select
    end subroutine reconstruction
    
    ! ===================================================================
    ! Flux Functions
    ! ===================================================================
    
    ! Rusanov flux
    subroutine rusanov_flux(q_face_left, q_face_right, flux, cfd)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdType), intent(in) :: cfd
        
        integer :: i
        real(dp) :: u_L, u_R, F_L, F_R, c_L, c_R, Smax
        
        c_L = cfd%config%wave_speed
        c_R = cfd%config%wave_speed
        
        do i = 1, cfd%domain%mesh%nnodes
            u_L = q_face_left(i)
            u_R = q_face_right(i)
            F_L = c_L * u_L
            F_R = c_R * u_R
            Smax = max(abs(c_L), abs(c_R))
            flux(i) = 0.5_dp * (F_L + F_R) - 0.5_dp * Smax * (u_R - u_L)
        end do
    end subroutine rusanov_flux
    
    ! Engquist-Osher flux
    subroutine engquist_osher_flux(q_face_left, q_face_right, flux, cfd)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdType), intent(in) :: cfd
        
        integer :: i
        real(dp) :: c, cp, cm, u_L, u_R
        
        c = cfd%config%wave_speed
        
        do i = 1, cfd%domain%mesh%nnodes
            cp = 0.5_dp * (c + abs(c))
            cm = 0.5_dp * (c - abs(c))
            u_L = q_face_left(i)
            u_R = q_face_right(i)
            flux(i) = cp * u_L + cm * u_R
        end do
    end subroutine engquist_osher_flux
    
    ! Inviscid flux selection
    subroutine inviscid_flux(q_face_left, q_face_right, flux, cfd)
        real(dp), intent(in) :: q_face_left(:), q_face_right(:)
        real(dp), intent(out) :: flux(:)
        type(CfdType), intent(in) :: cfd
        
        if (cfd%config%flux_type == 0) then
            call rusanov_flux(q_face_left, q_face_right, flux, cfd)
        else
            call engquist_osher_flux(q_face_left, q_face_right, flux, cfd)
        end if
    end subroutine inviscid_flux
    
    ! ===================================================================
    ! Residual Computation
    ! ===================================================================
    
    ! Compute residual (flux divergence) - CORRECTED VERSION
    subroutine residual(q, cfd)
        real(dp), intent(in) :: q(:)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i
        
        ! Apply boundary conditions first
        call boundary(cfd%solution%u, cfd)
        
        ! Reconstruction
        call reconstruction(q, cfd)
        
        ! Compute fluxes
        call inviscid_flux(cfd%solution%q_face_left, cfd%solution%q_face_right, &
                          cfd%solution%flux, cfd)
        
        ! Compute residual - corrected indexing
        do i = 1, cfd%domain%mesh%ncells
            cfd%solution%res(i) = -(cfd%solution%flux(i+1) - cfd%solution%flux(i)) / &
                                   cfd%domain%mesh%dx
        end do
    end subroutine residual
    
    ! ===================================================================
    ! Time Integration
    ! ===================================================================
    
    ! Update old field
    subroutine update_oldfield(qn, q, n)
        real(dp), intent(out) :: qn(:)
        real(dp), intent(in) :: q(:)
        integer, intent(in) :: n
        
        qn(1:n) = q(1:n)
    end subroutine update_oldfield
    
    ! 1st-order Runge-Kutta (Euler) - SIMPLIFIED
    subroutine runge_kutta_1(cfd)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i, j
        real(dp) :: dt
        
        dt = cfd%config%dt
        
        ! Apply boundary conditions
        call boundary(cfd%solution%u, cfd)
        
        ! Compute residual
        call residual(cfd%solution%u, cfd)
        
        ! Update solution
        do i = cfd%domain%ist, cfd%domain%ied - 1
            j = i - cfd%domain%ist + 1
            cfd%solution%u(i) = cfd%solution%u(i) + dt * cfd%solution%res(j)
        end do
        
        ! Apply boundary conditions again
        call boundary(cfd%solution%u, cfd)
        
        ! Save old solution
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine runge_kutta_1
    
    ! 2nd-order Runge-Kutta (Heun)
    subroutine runge_kutta_2(cfd)
        type(CfdType), intent(inout) :: cfd
        
        integer :: i, j
        real(dp) :: dt
        
        dt = cfd%config%dt
        
        ! Stage 1
        call boundary(cfd%solution%u, cfd)
        call residual(cfd%solution%u, cfd)
        
        do i = cfd%domain%ist, cfd%domain%ied - 1
            j = i - cfd%domain%ist + 1
            cfd%solution%u(i) = cfd%solution%u(i) + dt * cfd%solution%res(j)
        end do
        call boundary(cfd%solution%u, cfd)
        
        ! Stage 2
        call residual(cfd%solution%u, cfd)
        do i = cfd%domain%ist, cfd%domain%ied - 1
            j = i - cfd%domain%ist + 1
            cfd%solution%u(i) = 0.5_dp * cfd%solution%un(i) + &
                                0.5_dp * cfd%solution%u(i) + &
                                0.5_dp * dt * cfd%solution%res(j)
        end do
        call boundary(cfd%solution%u, cfd)
        
        call update_oldfield(cfd%solution%un, cfd%solution%u, cfd%domain%ntcells)
    end subroutine runge_kutta_2
    
    ! Runge-Kutta selection
    subroutine runge_kutta(cfd)
        type(CfdType), intent(inout) :: cfd
        
        select case(cfd%config%rk_order)
        case(1)
            call runge_kutta_1(cfd)
        case(2)
            call runge_kutta_2(cfd)
        case default
            call runge_kutta_1(cfd)
        end select
    end subroutine runge_kutta
    
    ! ===================================================================
    ! Simulation Driver
    ! ===================================================================
    
    ! Run simulation to final time
    function run_simulation(cfd, final_time) result(u_result)
        type(CfdType), intent(inout) :: cfd
        real(dp), intent(in) :: final_time
        real(dp), allocatable :: u_result(:)
        
        real(dp) :: t, dt, dt_old
        integer :: step, max_steps
        
        allocate(u_result(cfd%domain%mesh%ncells))
        
        t = 0.0_dp
        dt_old = cfd%config%dt
        dt = dt_old
        max_steps = 10000  ! Safety limit
        
        print *, "Starting time integration..."
        print *, "  Final time: ", final_time
        print *, "  Time step: ", dt
        print *, "  CFL number: ", cfd%config%cfl
        
        step = 0
        do while (t < final_time - 1.0e-12_dp .and. step < max_steps)
            step = step + 1
            
            if (t + dt > final_time) then
                dt = final_time - t
            end if
            
            cfd%config%dt = dt
            call runge_kutta(cfd)
            t = t + dt
            
            ! Progress report
            if (mod(step, 100) == 0) then
                print *, "  Step ", step, ", Time = ", t
            end if
        end do
        
        if (step >= max_steps) then
            print *, "Warning: Reached maximum number of steps (", max_steps, ")"
        end if
        
        cfd%config%dt = dt_old
        
        print *, "Time integration completed:"
        print *, "  Total steps: ", step
        print *, "  Final time: ", t
        
        ! Extract physical solution (without ghost cells)
        u_result = cfd%solution%u(cfd%domain%ist:cfd%domain%ied)
    end function run_simulation
    
    ! ===================================================================
    ! Main Analysis Function
    ! ===================================================================
    
    ! Perform ENO-WENO comparative analysis
    subroutine performEnoWenoAnalysis()
        type(CfdConfigType) :: config_eno3, config_weno3
        type(MeshType) :: mesh
        type(ComputationalDomainType) :: domain_eno3, domain_weno3
        type(CfdType) :: cfd_eno3, cfd_weno3
        real(dp), allocatable :: u_eno(:), u_weno(:), u_analytical(:)
        real(dp), allocatable :: xcc(:)
        integer :: i, ncells, iunit
        
        ! Initialize mesh
        call mesh%init()
        ncells = mesh%ncells
        allocate(xcc(ncells))
        xcc = mesh%xcc
        
        print *, "=========================================="
        print *, "Mesh parameters:"
        print *, "  ncells = ", ncells
        print *, "  dx = ", mesh%dx
        print *, "  L = ", mesh%L
        print *, "=========================================="
        
        ! Configure ENO3 - using simple 2nd order for stability
        config_eno3%recon_scheme = "eno"
        config_eno3%spatial_order = 3
        config_eno3%flux_type = 0
        config_eno3%rk_order = 1
        config_eno3%wave_speed = 1.0_dp
        config_eno3%final_time = 0.625_dp
        config_eno3%cfl = 1.0_dp
        config_eno3%dt = 0.0025_dp
        
        ! Configure WENO3
        config_weno3%recon_scheme = "weno"
        config_weno3%spatial_order = 3
        config_weno3%flux_type = 0
        config_weno3%rk_order = 1
        config_weno3%wave_speed = 1.0_dp
        config_weno3%final_time = 0.625_dp
        config_weno3%cfl = 1.0_dp 
        config_weno3%dt = 0.0025_dp
        
        ! Create domains
        call domain_eno3%init(mesh, config_eno3)
        call domain_weno3%init(mesh, config_weno3)
        
        ! Create CFD solvers
        call cfd_eno3%init(config_eno3, domain_eno3)
        call cfd_weno3%init(config_weno3, domain_weno3)
        
        ! Allocate arrays
        allocate(u_eno(ncells), u_weno(ncells), u_analytical(ncells))
        
        ! Run ENO simulation
        print *, "=========================================="
        print *, "Running ENO simulation..."
        print *, "  Scheme: ENO", config_eno3%spatial_order
        print *, "  Time step: ", config_eno3%dt
        print *, "=========================================="
        
        call init_field(cfd_eno3)
        u_eno = run_simulation(cfd_eno3, config_eno3%final_time)
        
        ! Run WENO simulation
        print *, "=========================================="
        print *, "Running WENO simulation..."
        print *, "  Scheme: WENO", config_weno3%spatial_order
        print *, "  Time step: ", config_weno3%dt
        print *, "=========================================="
        
        call init_field(cfd_weno3)
        u_weno = run_simulation(cfd_weno3, config_weno3%final_time)
        
        ! Compute analytical solution
        print *, "Computing analytical solution..."
        do i = 1, ncells
            u_analytical(i) = analytical_solution(xcc(i), config_weno3%final_time, &
                                                config_weno3%wave_speed, mesh%L)
        end do
        
        ! Write results to files
        print *, "Writing results to files..."
        
        ! Write ENO results
        open(newunit=iunit, file='eno_results.txt', status='replace')
        write(iunit, '(A)') '# x, u (ENO)'
        do i = 1, ncells
            write(iunit, '(2F12.6)') xcc(i), u_eno(i)
        end do
        close(iunit)
        
        ! Write WENO results
        open(newunit=iunit, file='weno_results.txt', status='replace')
        write(iunit, '(A)') '# x, u (WENO)'
        do i = 1, ncells
            write(iunit, '(2F12.6)') xcc(i), u_weno(i)
        end do
        close(iunit)
        
        ! Write analytical results
        open(newunit=iunit, file='analytical_results.txt', status='replace')
        write(iunit, '(A)') '# x, u (Analytical)'
        do i = 1, ncells
            write(iunit, '(2F12.6)') xcc(i), u_analytical(i)
        end do
        close(iunit)
        
        ! Print some statistics
        print *, "=========================================="
        print *, "Simulation statistics:"
        print *, "  ENO min/max: ", minval(u_eno), maxval(u_eno)
        print *, "  WENO min/max: ", minval(u_weno), maxval(u_weno)
        print *, "  Analytical min/max: ", minval(u_analytical), maxval(u_analytical)
        print *, "=========================================="
        
        print *, "Simulation completed successfully!"
        print *, "Results written to:"
        print *, "  eno_results.txt"
        print *, "  weno_results.txt"
        print *, "  analytical_results.txt"
        print *, ""
        print *, "To generate the comparison plot, run:"
        print *, "  python postprocess.py"
        print *, "=========================================="
        
        deallocate(u_eno, u_weno, u_analytical, xcc)
    end subroutine performEnoWenoAnalysis

end module cfd_solver

! ===================================================================
! Main Program
! ===================================================================
program main
    use cfd_solver
    implicit none
    
    print *, "=========================================="
    print *, "OneFLOW-CFD Solver for 1D Convection"
    print *, "ENO vs WENO Comparison"
    print *, "=========================================="
    
    call performEnoWenoAnalysis()
    
    print *, "Program finished successfully!"
    
end program main