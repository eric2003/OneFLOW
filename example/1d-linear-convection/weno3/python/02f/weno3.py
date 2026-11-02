import numpy as np
import matplotlib.pyplot as plt
import inflect

def initial_condition(x):
    """Initial condition: step function from 1.0 to 2.0 in [0.5, 1.0]"""
    u0 = np.zeros_like(x)
    for i in range(len(x)):
        if 0.5 <= x[i] <= 1.0:
            u0[i] = 2.0
        else:
            u0[i] = 1.0
    return u0

def analytical_solution(x, t, a, L):
    """Analytical solution with periodic boundary conditions"""
    x_shifted = (x - a * t + L) % L
    return initial_condition(x_shifted)    

# Compute residual (flux divergence) for all cells
def residual(q, cfd):
    reconstruction(q, cfd)
    inviscid_flux(cfd.q_face_left, cfd.q_face_right, cfd.flux, cfd)
    for i in range(cfd.ncells):
        cfd.res[i] = -(cfd.flux[i+1] - cfd.flux[i]) / cfd.mesh.dx
        
# Choose reconstruction method based on solver setting
def reconstruction(q, cfd):
    if cfd.config.reconstruction_scheme == 0:
        EnoReconstruction(q, cfd)
    elif cfd.config.reconstruction_scheme == 1:
        WenoReconstruction(q, cfd)

# --------------------------------------------------------------------------- #
# Reconstruction methods
# --------------------------------------------------------------------------- #          
def EnoReconstruction(q, cfd):
    """ENO reconstruction of interface values"""
    # Choose stencil by ENO method based on smoothest polynomial
    cfd.dd[0, :] = q
    
    # Compute divided differences
    for m in range(1, cfd.spatial_order):
        for j in range(0, cfd.ntcells-1):
            cfd.dd[m, j] = cfd.dd[m-1, j+1] - cfd.dd[m-1, j]
            
    # Select left-biased stencil for each node
    for i in range(cfd.nnodes):
        cfd.il[i] = i - 1
        for m in range(1, cfd.spatial_order):
            if abs(cfd.dd[m, cfd.il[i]-1+cfd.ishift]) <= abs(cfd.dd[m, cfd.il[i]+cfd.ishift]):
                cfd.il[i] -= 1
                
    # Select right-biased stencil for each node
    for i in range(cfd.nnodes):
        cfd.ir[i] = i
        for m in range(1, cfd.spatial_order):
            if abs(cfd.dd[m, cfd.ir[i]-1+cfd.ishift]) <= abs(cfd.dd[m, cfd.ir[i]+cfd.ishift]):
                cfd.ir[i] -= 1
    
    # Reconstruct values at cell interfaces (j+1/2)
    for i in range(cfd.nnodes):
        k1 = cfd.il[i]
        k2 = cfd.ir[i]
        l1 = i - k1
        l2 = i - k2
        cfd.q_face_left[i] = 0
        cfd.q_face_right[i] = 0
        for m in range(cfd.spatial_order):
            cfd.q_face_left[i] += q[k1 + cfd.ishift + m] * cfd.coef[l1, m]
            cfd.q_face_right[i] += q[k2 + cfd.ishift + m] * cfd.coef[l2, m]
            
# --------------------------------------------------------------------------- #
# Reconstruction methods
# --------------------------------------------------------------------------- #          
def EnoReconstructionOld(q, cfd):
    """ENO reconstruction of interface values"""
    # Choose stencil by ENO method based on smoothest polynomial
    cfd.dd[0, :] = q
    
    # Compute divided differences
    for m in range(1, cfd.spatial_order):
        for j in range(cfd.ntcells-m):
            cfd.dd[m, j] = cfd.dd[m-1, j+1] - cfd.dd[m-1, j]
            
    # Select left-biased stencil for each node
    for i in range(cfd.ist-1,cfd.ied+1):
        cfd.lmc[i] = i
        for m in range(1, cfd.spatial_order):
            if abs(cfd.dd[m, cfd.lmc[i]-1]) < abs(cfd.dd[m, cfd.lmc[i]]):
                cfd.lmc[i] -= 1
                
    # Reconstruct values at cell interfaces (j+1/2)
    for i in range(cfd.ist,cfd.ied+1):
        j = i - cfd.ist
        k1 = cfd.lmc[i-1]
        k2 = cfd.lmc[i  ]
        r1 = i-1 - k1
        r2 = i   - k2
        print(f"i,k1,k2,r1,r2={i,k1,k2,r1,r2}")
        cfd.q_face_left[j] = 0
        cfd.q_face_right[j] = 0
        for m in range(cfd.spatial_order):
            cfd.q_face_left[j] += q[k1 + m] * cfd.coef[r1+1, m]
            cfd.q_face_right[j] += q[k2 + m] * cfd.coef[r2, m]
            
def wc3L(v1,v2,v3):
    """WENO-3 nonlinear weights for left-biased stencil"""
    eps = 1.0e-6

    # Smoothness indicators
    s0 = (v3-v2)**2
    s1 = (v2-v1)**2

    # Compute nonlinear weights w0, w1
    d0 = 2.0/3.0
    d1 = 1.0/3.0
    
    c0 = d0 / ( (eps+s0)**2 )
    c1 = d1 / ( (eps+s1)**2 )

    w0 = c0 / ( c0 + c1 )
    w1 = c1 / ( c0 + c1 )

    # Candidate stencils
    q0 =  0.5 * v2 + 0.5 * v3
    q1 = -0.5 * v1 + 1.5 * v2

    # Reconstructed value at interface
    f = ( w0*q0 + w1*q1 )

    return f

def wc3R(v1,v2,v3):
    """WENO-3 nonlinear weights for right-biased stencil"""
    eps = 1.0e-6

    # Smoothness indicators
    s0 = (v2-v1)**2
    s1 = (v3-v2)**2

    # Compute nonlinear weights w0, w1
    d0 = 2.0/3.0
    d1 = 1.0/3.0
    
    c0 = d0 / ( (eps+s0)**2 )
    c1 = d1 / ( (eps+s1)**2 )

    w0 = c0 / ( c0 + c1 )
    w1 = c1 / ( c0 + c1 )

    # Candidate stencils
    q0 = 0.5 * v1 + 0.5 * v2
    q1 = 1.5 * v2 - 0.5 * v3

    # Reconstructed value at interface
    f = ( w0*q0 + w1*q1 )

    return f
    
# 3rd-order WENO reconstruction for left interface with periodic boundary
def weno3L_periodic(cfd,u,f):
    # i: ist-1, ist, ..., ied
    # j: 0, 1, ..., nx
    for i in range(cfd.ist - 1, cfd.ied):
        j = i - cfd.ist + 1
        v1 = u[i-1]
        v2 = u[i  ]
        v3 = u[i+1]
        f[j] = wc3L(v1,v2,v3)
        
# 3rd-order WENO reconstruction for right interface with periodic boundary
def weno3R_periodic(cfd,u,f):
    # i: ist, ist+1, ..., ied, ied+1
    # j: 0, 1, ..., nx
    for i in range(cfd.ist, cfd.ied + 1):
        j = i - cfd.ist
        v1 = u[i-1]
        v2 = u[i  ]
        v3 = u[i+1]
        f[j] = wc3R(v1,v2,v3)
        
# WENO (Weighted Essentially Non-Oscillatory) reconstruction
def WenoReconstruction(q, cfd):
    # Reconstruct values at cell interfaces (j+1/2)
    weno3L_periodic( cfd, q, cfd.q_face_left )
    weno3R_periodic( cfd, q, cfd.q_face_right )
    
fluxnames = [
        'Rusanov',
        'Engquist-Osher',
]

# Compute inviscid flux using selected Riemann solver
def inviscid_flux(q_face_left, q_face_right, flux, cfd):
    if cfd.config.flux_type == 0:
        rusanov_flux(q_face_left, q_face_right, flux, cfd)
    else:
        engquist_osher_flux(q_face_left, q_face_right, flux, cfd)
        
# --------------------------------------------------------------------------- #
# Numerical fluxes
# --------------------------------------------------------------------------- #
def rusanov_flux(q_face_left, q_face_right, flux, cfd):
    """Rusanov (local Lax-Friedrichs) flux"""
    for i in range(cfd.nnodes):
        u_L = q_face_left[i]
        u_R = q_face_right[i]
        F_L = cfd.config.wave_speed * u_L  # Flux from left state
        F_R = cfd.config.wave_speed * u_R  # Flux from right state
        alpha = abs(cfd.config.wave_speed) # Maximum wave speed
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)
        
def engquist_osher_flux(q_face_left, q_face_right, flux, cfd):
    """Engquist-Osher flux for linear convection"""
    for i in range(cfd.nnodes):
        c = cfd.config.wave_speed
        cp = 0.5*(c + abs(c))
        cm = 0.5*(c - abs(c))
        u_L = q_face_left[i]
        u_R = q_face_right[i]
       
        flux[i] = cp * u_L + cm * u_R

def periodic_boundary(u, cfd):
    """Apply periodic boundary conditions"""
    # Left ghost cells = right interior cells
    for ig in range(cfd.ighost):
        u[cfd.ist - 1 - ig] = u[cfd.ied - 1 - ig]
    
    # Right ghost cells = left interior cells
    for ig in range(cfd.ighost):
        u[cfd.ied + ig] = u[cfd.ist + ig]
    
    
# Apply periodic boundary conditions
def boundary(u, cfd):
    periodic_boundary(u, cfd)

# Copy current solution to old solution array
def update_oldfield(qn, q):
    qn[:] = q[:]

# Initialize reconstruction coefficients for different orders
def init_coef( spatial_order, coef ):
    if spatial_order == 1:
        coef[0] = [1.0]
        coef[1] = [1.0]
    elif spatial_order == 2:
        coef[0] = [3.0/2.0, -1.0/2.0]
        coef[1] = [1.0/2.0,  1.0/2.0]
        coef[2] = [-1.0/2.0, 3.0/2.0]
    elif spatial_order == 3:
        coef[0] = [ 11.0/6.0, -7.0/6.0,  1.0/3.0 ]
        coef[1] = [  1.0/3.0,  5.0/6.0, -1.0/6.0 ]
        coef[2] = [ -1.0/6.0,  5.0/6.0,  1.0/3.0 ]
        coef[3] = [  1.0/3.0, -7.0/6.0, 11.0/6.0 ]
    elif spatial_order == 4:
        coef[0] = [ 25.0/12.0, -23.0/12.0,  13.0/12.0,  -1.0/4.0 ]
        coef[1] = [   1.0/4.0,  13.0/12.0,  -5.0/12.0,  1.0/12.0 ]
        coef[2] = [ -1.0/12.0,   7.0/12.0,   7.0/12.0, -1.0/12.0 ]
        coef[3] = [  1.0/12.0,  -5.0/12.0,  13.0/12.0,   1.0/4.0 ]
        coef[4] = [  -1.0/4.0,  13.0/12.0, -23.0/12.0, 25.0/12.0 ]
    elif spatial_order == 5:
        coef[0] = [ 137.0/60.0, -163.0/60.0, 137.0/60.0,  -21.0/20.0,    1.0/5.0 ]
        coef[1] = [    1.0/5.0,   77.0/60.0, -43.0/60.0,   17.0/60.0,  -1.0/20.0 ]
        coef[2] = [  -1.0/20.0,    9.0/20.0,  47.0/60.0,  -13.0/60.0,   1.0/30.0 ]
        coef[3] = [   1.0/30.0,  -13.0/60.0,  47.0/60.0,    9.0/20.0,  -1.0/20.0 ]
        coef[4] = [  -1.0/20.0,   17.0/60.0, -43.0/60.0,   77.0/60.0,    1.0/5.0 ]
        coef[5] = [    1.0/5.0,  -21.0/20.0, 137.0/60.0, -163.0/60.0, 137.0/60.0 ]
    elif spatial_order == 6:
        coef[0] = [ 49.0/20.0, -71.0/20.0,   79.0/20.0, -163.0/60.0,  31.0/30.0,  -1.0/6.0 ]
        coef[1] = [   1.0/6.0,  29.0/20.0,  -21.0/20.0,   37.0/60.0, -13.0/60.0,  1.0/30.0 ]
        coef[2] = [ -1.0/30.0,  11.0/30.0,   19.0/20.0,  -23.0/60.0,   7.0/60.0, -1.0/60.0 ]
        coef[3] = [  1.0/60.0,  -2.0/15.0,   37.0/60.0,   37.0/60.0,  -2.0/15.0,  1.0/60.0 ]
        coef[4] = [ -1.0/60.0,   7.0/60.0,  -23.0/60.0,   19.0/20.0,  11.0/30.0, -1.0/30.0 ]
        coef[5] = [  1.0/30.0, -13.0/60.0,   37.0/60.0,  -21.0/20.0,  29.0/20.0,   1.0/6.0 ]
        coef[6] = [  -1.0/6.0,  31.0/30.0, -163.0/60.0,   79.0/20.0, -71.0/20.0, 49.0/20.0 ]
    elif spatial_order == 7:
        coef[0] = [ 363.0/140.0, -617.0/140.0,  853.0/140.0, -2341.0/420.0,  667.0/210.0,   -43.0/42.0,     1.0/7.0 ]
        coef[1] = [     1.0/7.0,  223.0/140.0, -197.0/140.0,   153.0/140.0, -241.0/420.0,   37.0/210.0,   -1.0/42.0 ]
        coef[2] = [   -1.0/42.0,    13.0/42.0,  153.0/140.0,  -241.0/420.0,  109.0/420.0,  -31.0/420.0,   1.0/105.0 ]
        coef[3] = [   1.0/105.0,  -19.0/210.0,  107.0/210.0,   319.0/420.0, -101.0/420.0,     5.0/84.0,  -1.0/140.0 ]
        coef[4] = [  -1.0/140.0,     5.0/84.0, -101.0/420.0,   319.0/420.0,  107.0/210.0,  -19.0/210.0,   1.0/105.0 ]
        coef[5] = [   1.0/105.0,  -31.0/420.0,  109.0/420.0,  -241.0/420.0,  153.0/140.0,    13.0/42.0,   -1.0/42.0 ]
        coef[6] = [   -1.0/42.0,   37.0/210.0, -241.0/420.0,   153.0/140.0, -197.0/140.0,  223.0/140.0,     1.0/7.0 ]
        coef[7] = [     1.0/7.0,   -43.0/42.0,  667.0/210.0, -2341.0/420.0,  853.0/140.0, -617.0/140.0, 363.0/140.0 ]

# Initialize flow field with piecewise constant distribution
def init_field(cfd):
    for i in range(cfd.ist, cfd.ied):
        j = i - cfd.ist
        if 0.5 <= cfd.mesh.xcc[j] <= 1.0:
            cfd.u[i] = 2.0
        else:
            cfd.u[i] = 1.0
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)
    
# Select Runge-Kutta time integration scheme
def runge_kutta(cfd):
    rk_order = cfd.config.rk_order
    if rk_order == 1:
        runge_kutta_1(cfd)
    elif rk_order == 2:
        runge_kutta_2(cfd)
    else:
        runge_kutta_3(cfd)
    
# 1st-order explicit Euler time integration
def runge_kutta_1(cfd):
    dt = cfd.config.dt
    residual(cfd.u, cfd)
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)
    
# 2nd-order Runge-Kutta (Heun's method) time integration
def runge_kutta_2(cfd):
    dt = cfd.config.dt
    residual(cfd.u, cfd)
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    
    residual(cfd.u, cfd)
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = 0.5 * cfd.un[j] + 0.5 * cfd.u[j] + 0.5 * dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)

# 3rd-order Runge-Kutta (SSPRK3) time integration
def runge_kutta_3(cfd):
    dt = cfd.config.dt
    residual(cfd.u, cfd)
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    
    residual(cfd.u, cfd)
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = 0.75 * cfd.un[j] + 0.25 * cfd.u[j] + 0.25 * dt * cfd.res[i]
    boundary(cfd.u, cfd)
    
    residual(cfd.u, cfd)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = c1 * cfd.un[j] + c2 * cfd.u[j] + c3 * dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)

# Visualize numerical and analytical solutions
def visualize(cfd):
    with open('solution.plt', 'w') as f:
        for i in range(cfd.ist, cfd.ied):
            j = i - cfd.ist
            f.write(f"{cfd.mesh.xcc[j]:20.10e}{cfd.u[i]:20.10e}\n")
            
    # Visualization setup
    u_numerical = np.copy(cfd.u[cfd.ist:cfd.ied])
    print(f'u_numerical.size={u_numerical.size}')
    print(f'cfd.mesh.xcc.size={cfd.mesh.xcc.size}')
    # Compute analytical solution
    u_analytical = analytical_solution(cfd.mesh.xcc, cfd.config.final_time, cfd.config.wave_speed, cfd.mesh.L)
    print(f'u_analytical.size={u_analytical.size}')
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.scatter(cfd.mesh.xcc, u_numerical, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=f'Numerical (Rusanov)')
    plt.plot(cfd.mesh.xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    
    p = inflect.engine()
    iorder_str = p.ordinal(cfd.spatial_order)
    rk_str = p.ordinal(cfd.config.rk_order)
    plt.title(f'1D Convection Equation at t = {cfd.config.final_time:.3f} using {iorder_str}-order WENO and {rk_str}-order Runge-Kutta methods')
    plt.legend()
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)    
    plt.tight_layout()
    plt.show()

# Mesh class: defines computational grid
class Mesh:
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 2.0
        self.ncells = 40
        self.nnodes = self.ncells + 1
        self.nx = self.ncells
        self.x = np.zeros(self.nnodes)
        self.xcc = np.zeros(self.ncells)
        self.init_mesh()
        
    def init_mesh(self):
        self.L = self.xmax - self.xmin
        self.dx = self.L / self.ncells
        
        # Generate node coordinates
        for i in range(self.nnodes):
            self.x[i] = self.xmin + i * self.dx
            
        # Generate cell center coordinates
        for i in range(self.ncells):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i+1])

class SimulationConfig:
    def __init__(self):
        self.reconstruction_scheme = 0  # 0=ENO, 1=WENO
        self.flux_type = 0      # 0=Rusanov, 1=Engquist-Osher
        self.rk_order = 1
        self.wave_speed = 1.0
        self.final_time = 0.625
        self.dt = 0.025
        
        
# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, config, mesh, spatial_order):
        self.config = config
        self.mesh = mesh
        self.nx = mesh.nx
        self.spatial_order = spatial_order
        self.ighost = spatial_order  # Number of ghost cells
        self.ishift = self.ighost + 1
        #self.ishift = self.ighost
        self.ncells = mesh.ncells
        self.nnodes = mesh.nnodes
        self.ist = 0 + self.ishift  # Start index of physical cells
        self.ied = self.ncells + self.ist  # End index of physical cells
        self.ntcells = self.ncells + 2 * self.ishift  # Total cells including ghost regions
        print(f"self.ncells={self.ncells}")
        print(f"self.spatial_order={self.spatial_order}")
        print(f"self.ighost={self.ighost}")
        print(f"self.ishift={self.ishift}")
        print(f"self.ist={self.ist}")
        print(f"self.ied={self.ied}")

        # Stencil selection arrays
        self.il = np.zeros(self.nnodes, dtype=int)
        self.ir = np.zeros(self.nnodes, dtype=int)
        self.lmc = np.zeros(self.ntcells, dtype=int)
        
        # Reconstruction coefficients and divided differences
        self.coef = np.zeros((spatial_order + 1, spatial_order))
        self.dd = np.zeros((spatial_order, self.ntcells))
        
        # Interface values and fluxes
        self.q_face_left = np.zeros(self.nnodes)  # Left interface value
        self.q_face_right = np.zeros(self.nnodes)  # Right interface value
        self.flux = np.zeros(self.nnodes)
        self.res = np.zeros(self.ncells)  # Residual array

        # Solution arrays
        self.u = np.zeros(self.ntcells)  # Current solution
        self.un = np.zeros(self.ntcells)  # Previous time step solution
        
        init_coef(self.spatial_order, self.coef)
        
# --------------------------------------------------------------------------- #
# Simulation runners
# --------------------------------------------------------------------------- #
def run_simulation(cfd, final_time):
    t = 0.0
    dt_old = cfd.config.dt
    dt = dt_old
    while t < final_time:
        if t + dt > final_time:
            dt = final_time - t
        cfd.config.dt = dt  # temporary adjustment for last step
        runge_kutta(cfd)
        t += dt
    cfd.config.dt = dt_old
    return cfd.u[cfd.ist:cfd.ied].copy()        

# Perform ENO-WENO comparative analysis
def performEnoWenoAnalysis():
    mesh = Mesh()
    config = SimulationConfig()
    u_analytical = analytical_solution(mesh.xcc, config.final_time, config.wave_speed, mesh.L)
    
    config.rk_order = 1
    config.dt = 0.0025
    
    u_list = []
    # ENO
    config.reconstruction_scheme = 0
    cfd = Cfd(config, mesh, spatial_order=3)
    init_field(cfd)
    u_eno = run_simulation(cfd, config.final_time)
    u_list.append(u_eno)
    
    # WENO
    config.reconstruction_scheme = 1
    cfd = Cfd(config, mesh, spatial_order=3)
    init_field(cfd)
    u_weno = run_simulation(cfd, config.final_time)
    u_list.append(u_weno)

    plot_EnoWeno_Analysis(config, mesh.xcc, u_list, u_analytical)
    
# Plot ENO-WENO comparison results
def plot_EnoWeno_Analysis(config, xcc, u_list, u_analytical):
    # Define line styles with different colors and markers
    styles = [
        {'color': 'black', 'linestyle': '-', 'marker': 'o'},
        {'color': 'blue', 'linestyle': '--', 'marker': 's'},
        {'color': 'black', 'linestyle': '-', 'marker': '^'},
        {'color': 'blue', 'linestyle': '--', 'marker': 'v'},
        {'color': 'black', 'linestyle': '-', 'marker': '<'},
        {'color': 'blue', 'linestyle': '--', 'marker': '>'},
        {'color': 'black', 'linestyle': '-', 'marker': 'D'},
    ]

    n = len(u_list)
    num_styles = len(styles)

    p = inflect.engine()
    rk_str = p.ordinal(config.rk_order)    
    
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.title(f'1D Convection Equation at t = {config.final_time:.3f} using 3rd-order ENO&WENO and {rk_str}-order Runge-Kutta methods')
    for i in range(0, n):
        if i == 0:
            lable = 'Numerical (Rusanov)ENO3'
        else:
            lable = 'Numerical (Rusanov)WENO3'
        style = styles[i % num_styles]
        plt.plot(xcc, u_list[i], marker=style['marker'], markerfacecolor='none', linestyle=style['linestyle'], color=style['color'], \
                markersize=5, linewidth=0.5, alpha=1.0, label=f'{lable}')
    plt.plot(xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    plt.show()    
    
# Main execution function
def main():
    performEnoWenoAnalysis()

if __name__ == "__main__":
    main()