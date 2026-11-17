import numpy as np
import matplotlib.pyplot as plt

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
    inviscid_flux(cfd.up1_2m, cfd.up1_2p, cfd.flux, cfd)
    for i in range(cfd.ncells):
        cfd.res[i] = -(cfd.flux[i+1] - cfd.flux[i]) / cfd.mesh.dx
        
# Choose reconstruction method based on solver setting
def reconstruction(q, cfd):
    if cfd.solver.interpolation == 0:
        EnoReconstruction(q, cfd)
    elif cfd.solver.interpolation == 1:
        WenoReconstruction(q, cfd)

# --------------------------------------------------------------------------- #
# Reconstruction methods
# --------------------------------------------------------------------------- #          
def EnoReconstruction(q, cfd):
    """ENO reconstruction of interface values"""
    # Choose stencil by ENO method based on smoothest polynomial
    cfd.dd[0, :] = q
    
    # Compute divided differences
    for m in range(1, cfd.iorder):
        for j in range(0, cfd.ntcells-1):
            cfd.dd[m, j] = cfd.dd[m-1, j+1] - cfd.dd[m-1, j]
            
    # Select left-biased stencil for each node
    for i in range(cfd.nnodes):
        cfd.il[i] = i - 1
        for m in range(1, cfd.iorder):
            if abs(cfd.dd[m, cfd.il[i]-1+cfd.ishift]) <= abs(cfd.dd[m, cfd.il[i]+cfd.ishift]):
                cfd.il[i] -= 1
                
    # Select right-biased stencil for each node
    for i in range(cfd.nnodes):
        cfd.ir[i] = i
        for m in range(1, cfd.iorder):
            if abs(cfd.dd[m, cfd.ir[i]-1+cfd.ishift]) <= abs(cfd.dd[m, cfd.ir[i]+cfd.ishift]):
                cfd.ir[i] -= 1
    
    # Reconstruct values at cell interfaces (j+1/2)
    for i in range(cfd.nnodes):
        k1 = cfd.il[i]
        k2 = cfd.ir[i]
        l1 = i - k1
        l2 = i - k2
        cfd.up1_2m[i] = 0
        cfd.up1_2p[i] = 0
        for m in range(cfd.iorder):
            cfd.up1_2m[i] += q[k1 + cfd.ishift + m] * cfd.coef[l1, m]
            cfd.up1_2p[i] += q[k2 + cfd.ishift + m] * cfd.coef[l2, m]
            
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
    weno3L_periodic( cfd, q, cfd.up1_2m )
    weno3R_periodic( cfd, q, cfd.up1_2p )
    
fluxnames = [
        'Rusanov',
        'Engquist-Osher',
]

# Compute inviscid flux using selected Riemann solver
def inviscid_flux(up1_2m, up1_2p, flux, cfd):
    if cfd.solver.iflux == 0:
        rusanov_flux(up1_2m, up1_2p, flux, cfd)
    else:
        engquist_osher_flux(up1_2m, up1_2p, flux, cfd)
        
# --------------------------------------------------------------------------- #
# Numerical fluxes
# --------------------------------------------------------------------------- #
def rusanov_flux(up1_2m, up1_2p, flux, cfd):
    """Rusanov (local Lax-Friedrichs) flux"""
    for i in range(cfd.nnodes):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        F_L = cfd.solver.c * u_L  # Flux from left state
        F_R = cfd.solver.c * u_R  # Flux from right state
        alpha = abs(cfd.solver.c) # Maximum wave speed
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)        
        
def engquist_osher_flux(up1_2m, up1_2p, flux, cfd):
    """Engquist-Osher flux for linear convection"""
    for i in range(cfd.nnodes):
        c = cfd.solver.c
        cp = 0.5*(c + abs(c))
        cm = 0.5*(c - abs(c))
        
        u_L = up1_2m[i]
        u_R = up1_2p[i]
       
       
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
def init_coef( iorder, coef ):
    if iorder == 1:
        coef[0] = [1.0]
        coef[1] = [1.0]
    elif iorder == 2:
        coef[0] = [3.0/2.0, -1.0/2.0]
        coef[1] = [1.0/2.0,  1.0/2.0]
        coef[2] = [-1.0/2.0, 3.0/2.0]
    elif iorder == 3:
        coef[0] = [ 11.0/6.0, -7.0/6.0,  1.0/3.0 ]
        coef[1] = [  1.0/3.0,  5.0/6.0, -1.0/6.0 ]
        coef[2] = [ -1.0/6.0,  5.0/6.0,  1.0/3.0 ]
        coef[3] = [  1.0/3.0, -7.0/6.0, 11.0/6.0 ]
    elif iorder == 4:
        coef[0] = [ 25.0/12.0, -23.0/12.0,  13.0/12.0,  -1.0/4.0 ]
        coef[1] = [   1.0/4.0,  13.0/12.0,  -5.0/12.0,  1.0/12.0 ]
        coef[2] = [ -1.0/12.0,   7.0/12.0,   7.0/12.0, -1.0/12.0 ]
        coef[3] = [  1.0/12.0,  -5.0/12.0,  13.0/12.0,   1.0/4.0 ]
        coef[4] = [  -1.0/4.0,  13.0/12.0, -23.0/12.0, 25.0/12.0 ]
    elif iorder == 5:
        coef[0] = [ 137.0/60.0, -163.0/60.0, 137.0/60.0,  -21.0/20.0,    1.0/5.0 ]
        coef[1] = [    1.0/5.0,   77.0/60.0, -43.0/60.0,   17.0/60.0,  -1.0/20.0 ]
        coef[2] = [  -1.0/20.0,    9.0/20.0,  47.0/60.0,  -13.0/60.0,   1.0/30.0 ]
        coef[3] = [   1.0/30.0,  -13.0/60.0,  47.0/60.0,    9.0/20.0,  -1.0/20.0 ]
        coef[4] = [  -1.0/20.0,   17.0/60.0, -43.0/60.0,   77.0/60.0,    1.0/5.0 ]
        coef[5] = [    1.0/5.0,  -21.0/20.0, 137.0/60.0, -163.0/60.0, 137.0/60.0 ]
    elif iorder == 6:
        coef[0] = [ 49.0/20.0, -71.0/20.0,   79.0/20.0, -163.0/60.0,  31.0/30.0,  -1.0/6.0 ]
        coef[1] = [   1.0/6.0,  29.0/20.0,  -21.0/20.0,   37.0/60.0, -13.0/60.0,  1.0/30.0 ]
        coef[2] = [ -1.0/30.0,  11.0/30.0,   19.0/20.0,  -23.0/60.0,   7.0/60.0, -1.0/60.0 ]
        coef[3] = [  1.0/60.0,  -2.0/15.0,   37.0/60.0,   37.0/60.0,  -2.0/15.0,  1.0/60.0 ]
        coef[4] = [ -1.0/60.0,   7.0/60.0,  -23.0/60.0,   19.0/20.0,  11.0/30.0, -1.0/30.0 ]
        coef[5] = [  1.0/30.0, -13.0/60.0,   37.0/60.0,  -21.0/20.0,  29.0/20.0,   1.0/6.0 ]
        coef[6] = [  -1.0/6.0,  31.0/30.0, -163.0/60.0,   79.0/20.0, -71.0/20.0, 49.0/20.0 ]
    elif iorder == 7:
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
    rk = cfd.solver.rk
    if rk == 1:
        runge_kutta_1(cfd)
    elif rk == 2:
        runge_kutta_2(cfd)
    else:
        runge_kutta_3(cfd)
    
# 1st-order explicit Euler time integration
def runge_kutta_1(cfd):
    dt = cfd.solver.dt
    residual(cfd.u, cfd)
    for i in range(cfd.ncells):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)
    
# 2nd-order Runge-Kutta (Heun's method) time integration
def runge_kutta_2(cfd):
    dt = cfd.solver.dt
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
    dt = cfd.solver.dt
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

# Get ordinal suffix for number formatting (st, nd, rd, th)
def get_ordinal_numbers(order):
    if order == 1:
        return 'st'
    elif order == 2:
        return 'nd'
    elif order == 3:
        return 'rd'
    else:
        return 'th'

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
    u_analytical = analytical_solution(cfd.mesh.xcc, cfd.solver.T, cfd.solver.c, cfd.mesh.L)
    print(f'u_analytical.size={u_analytical.size}')
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.scatter(cfd.mesh.xcc, u_numerical, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=f'Numerical (Rusanov)')
    plt.plot(cfd.mesh.xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    ordinal1 = get_ordinal_numbers(cfd.iorder)
    ordinal2 = get_ordinal_numbers(cfd.solver.rk)
    plt.title(f'1D Convection Equation at t = {cfd.solver.T:.3f} using {cfd.iorder}{ordinal1}-order WENO and {cfd.solver.rk}{ordinal2}-order Runge-Kutta methods')
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

# Solver class: stores numerical method parameters
class Solver:
    def __init__(self):
        # interpolation: 0 for ENO, 1 for WENO
        self.interpolation = 0
        self.iflux = 0
        self.rk = 1
        self.c = 1.0
        self.T = 0.625
        #self.dt = .0025
        self.dt = .025
        
# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, solver, mesh, iorder):
        self.solver = solver
        self.mesh = mesh
        self.nx = mesh.nx
        self.iorder = iorder
        self.ighost = iorder  # Number of ghost cells
        self.ishift = self.ighost + 1
        #self.ishift = self.ighost
        self.ncells = mesh.ncells
        self.nnodes = mesh.nnodes
        self.ist = 0 + self.ishift  # Start index of physical cells
        self.ied = self.ncells + self.ist  # End index of physical cells
        self.ntcells = self.ncells + 2 * self.ishift  # Total cells including ghost regions
        print(f"self.ncells={self.ncells}")
        print(f"self.iorder={self.iorder}")
        print(f"self.ighost={self.ighost}")
        print(f"self.ishift={self.ishift}")
        print(f"self.ist={self.ist}")
        print(f"self.ied={self.ied}")

        # Stencil selection arrays
        self.il = np.zeros(self.nnodes, dtype=int)
        self.ir = np.zeros(self.nnodes, dtype=int)
        
        # Reconstruction coefficients and divided differences
        self.coef = np.zeros((iorder + 1, iorder))
        self.dd = np.zeros((iorder, self.ntcells))
        
        # Interface values and fluxes
        self.up1_2m = np.zeros(self.nnodes)  # Left interface value
        self.up1_2p = np.zeros(self.nnodes)  # Right interface value
        self.flux = np.zeros(self.nnodes)
        self.res = np.zeros(self.ncells)  # Residual array

        # Solution arrays
        self.u = np.zeros(self.ntcells)  # Current solution
        self.un = np.zeros(self.ntcells)  # Previous time step solution
        
        init_coef(self.iorder, self.coef)
        
# --------------------------------------------------------------------------- #
# Simulation runners
# --------------------------------------------------------------------------- #
def run_simulation(cfd, final_time):
    t = 0.0
    dt_old = cfd.solver.dt
    dt = dt_old
    while t < final_time:
        if t + dt > final_time:
            dt = final_time - t
        cfd.solver.dt = dt  # temporary adjustment for last step
        runge_kutta(cfd)
        t += dt
    cfd.solver.dt = dt_old
    return cfd.u[cfd.ist:cfd.ied].copy()        

# Perform ENO-WENO comparative analysis
def performEnoWenoAnalysis():
    mesh = Mesh()
    solver = Solver()
    u_analytical = analytical_solution(mesh.xcc, solver.T, solver.c, mesh.L)
    
    solver.rk = 1
    solver.dt = 0.0025
    solver.iorder = 3
    
    u_list = []
    # ENO
    solver.interpolation = 0
    cfd = Cfd(solver, mesh, iorder=3)
    init_field(cfd)
    u_eno = run_simulation(cfd, solver.T)
    u_list.append(u_eno)
    
    # WENO
    solver.interpolation = 1
    cfd = Cfd(solver, mesh, iorder=3)
    init_field(cfd)
    u_weno = run_simulation(cfd, solver.T)
    u_list.append(u_weno)    
    
       
    plot_EnoWeno_Analysis(solver, mesh.xcc, u_list, u_analytical)
    
# Plot ENO-WENO comparison results
def plot_EnoWeno_Analysis(solver, xcc, u_list, u_analytical):
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
    ordinal = get_ordinal_numbers(solver.rk)
    
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.title(f'1D Convection Equation at t = {solver.T:.3f} using 3rd-order ENO&WENO and {solver.rk}{ordinal}-order Runge-Kutta methods')
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