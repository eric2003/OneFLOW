import numpy as np
import matplotlib.pyplot as plt

# 初始条件
def initial_condition(x):
    u0 = np.zeros_like(x)
    for i in range(len(x)):
        if 0.5 <= x[i] <= 1.0:
            u0[i] = 2.0
        else:
            u0[i] = 1.0
    return u0

# 理论解
def analytical_solution(x, t, a, L):
    # 初始条件沿 x - at 平移
    x_shifted = x - a * t
    return initial_condition((x_shifted + L) % L)  # 周期边界条件  

def residual(q, cfd):
    reconstruction(q, cfd)
    inviscid_flux(cfd.up1_2m, cfd.up1_2p, cfd.flux, cfd)
    for i in range(cfd.nx):
        cfd.res[i] = -(cfd.flux[i + 1] - cfd.flux[i]) / cfd.mesh.dx
        
def reconstruction(q, cfd):
    if cfd.solver.interpolation == 0:
        EnoReconstruction(q, cfd)
    elif cfd.solver.interpolation == 1:
        WenoReconstruction(q, cfd)
    
def EnoReconstruction(q, cfd):
    # Choose the stencil by ENO method
    cfd.dd[0, 0:cfd.ntcell-1] = q[0:cfd.ntcell-1]
    
    for m in range(1, cfd.iorder):
        for j in range(0, cfd.ntcell-1):
            cfd.dd[m, j] = cfd.dd[m-1, j+1] - cfd.dd[m-1, j]
            
    for i in range(cfd.nx + 1):
        cfd.il[i] = i - 1
        for m in range(1, cfd.iorder):
            if abs(cfd.dd[m, cfd.il[i]-1+cfd.ishift]) <= abs(cfd.dd[m, cfd.il[i]+cfd.ishift]):
                cfd.il[i] -= 1
                
    for i in range(cfd.nx + 1):
        cfd.ir[i] = i
        for m in range(1, cfd.iorder):
            if abs(cfd.dd[m, cfd.ir[i]-1+cfd.ishift]) <= abs(cfd.dd[m, cfd.ir[i]+cfd.ishift]):
                cfd.ir[i] -= 1
    
    # Reconstruction u(j+1/2)
    for i in range(cfd.nx + 1):
        k1 = cfd.il[i]
        k2 = cfd.ir[i]
        l1 = i - k1
        l2 = i - k2
        cfd.up1_2m[i] = 0
        cfd.up1_2p[i] = 0
        for m in range(cfd.iorder):
            cfd.up1_2m[i] += q[k1 + cfd.ishift + m] * cfd.coef[l1, m]
            cfd.up1_2p[i] += q[k2 + cfd.ishift + m] * cfd.coef[l2, m]
            
#---------------------------------------------------------------------------#
#nonlinear weights for upwind direction
#---------------------------------------------------------------------------#
def wc3L(v1,v2,v3):
    eps = 1.0e-6

    # smoothness indicators
    s0 = (v3-v2)**2
    s1 = (v2-v1)**2

    # computing nonlinear weights w1,w2
    d0 = 2.0/3.0
    d1 = 1.0/3.0
    
    c0 = d0 / ( (eps+s0)**2 )
    c1 = d1 / ( (eps+s1)**2 )

    w0 = c0 / ( c0 + c1 )
    w1 = c1 / ( c0 + c1 )

    # candiate stencils
    q0 =  0.5 * v2 + 0.5 * v3
    q1 = -0.5 * v1 + 1.5 * v2

    # reconstructed value at interface
    f = ( w0*q0 + w1*q1 )

    return f

#---------------------------------------------------------------------------#
#nonlinear weights for downwind direction
#---------------------------------------------------------------------------#
def wc3R(v1,v2,v3):
    eps = 1.0e-6

    # smoothness indicators
    s0 = (v2-v1)**2
    s1 = (v3-v2)**2

    # computing nonlinear weights w1,w2
    d0 = 2.0/3.0
    d1 = 1.0/3.0
    
    c0 = d0 / ( (eps+s0)**2 )
    c1 = d1 / ( (eps+s1)**2 )

    w0 = c0 / ( c0 + c1 )
    w1 = c1 / ( c0 + c1 )

    # candiate stencils
    q0 = 0.5 * v1 + 0.5 * v2
    q1 = 1.5 * v2 - 0.5 * v3

    # reconstructed value at interface
    f = ( w0*q0 + w1*q1 )

    return f
    
            
def weno3L_periodic(cfd,u,f):
    #i:ist-1,ist,...,ied
    #j:0,1,...,nx
    for i in range(cfd.ist - 1, cfd.ied + 1):
        j = i - cfd.ist + 1
        v1 = u[i-1]
        v2 = u[i  ]
        v3 = u[i+1]
        f[j] = wc3L(v1,v2,v3)
        
def weno3R_periodic(cfd,u,f):
    #i:ist,ist+1,...,ied,ied+1
    #j:0,1,...,nx
    for i in range(cfd.ist, cfd.ied + 2):
        j = i - cfd.ist
        v1 = u[i-1]
        v2 = u[i  ]
        v3 = u[i+1]
        f[j] = wc3R(v1,v2,v3)
        
            
def WenoReconstruction(q, cfd):
    # Reconstruction u(j+1/2)
    weno3L_periodic( cfd, q, cfd.up1_2m )
    weno3R_periodic( cfd, q, cfd.up1_2p )

fluxnames = [
        'Rusanov',
        'Engquist-Osher',
]

def inviscid_flux(up1_2m, up1_2p, flux, cfd):
    if cfd.solver.iflux == 0:
        rusanov_flux(up1_2m, up1_2p, flux, cfd)
    else:
        engquist_osher_flux(up1_2m, up1_2p, flux, cfd)

def engquist_osher_flux(up1_2m, up1_2p, flux, cfd):
    for i in range(cfd.nx + 1):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        
        cp = 0.5 * ( cfd.solver.c + abs(cfd.solver.c) )
        cm = 0.5 * ( cfd.solver.c - abs(cfd.solver.c) )
        
        flux[i] = cp * u_L + cm * u_R
        
def rusanov_flux(up1_2m, up1_2p, flux, cfd):
    for i in range(cfd.nx + 1):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        F_L = cfd.solver.c * u_L  # 左状态通量
        F_R = cfd.solver.c * u_R  # 右状态通量
        alpha = abs(cfd.solver.c) # 最大波速
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)
                

def boundary(u, cfd):
    for i in range(-cfd.ighost, 1):
        u[cfd.ist - 1 + i] = u[cfd.ied + i]
    for i in range(1, cfd.ighost + 2):
        u[cfd.ied + i] = u[cfd.ist - 1 + i]

def update_oldfield(qn, q):
    qn[:] = q[:]

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

def init_field(cfd):
    for i in range(cfd.ist, cfd.ied + 1):
        j = i - cfd.ist
        if 0.5 <= cfd.mesh.xcc[j] <= 1.0:
            cfd.u[i] = 2.0
        else:
            cfd.u[i] = 1.0
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)
    
def runge_kutta(cfd):
    rk = cfd.solver.rk
    if rk == 1:
        runge_kutta_1(cfd)
    elif rk == 2:
        runge_kutta_2(cfd)
    else:
        runge_kutta_3(cfd)
    
def runge_kutta_1(cfd):
    dt = cfd.solver.dt
    residual(cfd.u, cfd)
    for i in range(cfd.nx):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)
    
def runge_kutta_2(cfd):
    dt = cfd.solver.dt
    residual(cfd.u, cfd)
    for i in range(cfd.nx):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    
    residual(cfd.u, cfd)
    for i in range(cfd.nx):
        j = i + cfd.ishift
        cfd.u[j] = 0.5 * cfd.un[j] + 0.5 * cfd.u[j] + 0.5 * dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)

def runge_kutta_3(cfd):
    dt = cfd.solver.dt
    residual(cfd.u, cfd)
    for i in range(cfd.nx):
        j = i + cfd.ishift
        cfd.u[j] = cfd.u[j] + dt * cfd.res[i]
    boundary(cfd.u, cfd)
    
    residual(cfd.u, cfd)
    for i in range(cfd.nx):
        j = i + cfd.ishift
        cfd.u[j] = 0.75 * cfd.un[j] + 0.25 * cfd.u[j] + 0.25 * dt * cfd.res[i]
    boundary(cfd.u, cfd)
    
    residual(cfd.u, cfd)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in range(cfd.nx):
        j = i + cfd.ishift
        cfd.u[j] = c1 * cfd.un[j] + c2 * cfd.u[j] + c3 * dt * cfd.res[i]
    boundary(cfd.u, cfd)
    update_oldfield(cfd.un, cfd.u)

def get_ordinal_numbers(order):
    if order == 1:
        return 'st'
    elif order == 2:
        return 'nd'
    elif order == 3:
        return 'rd'
    else:
        return 'th'

def visualize(cfd):
    with open('solution.plt', 'w') as f:
        for i in range(cfd.ist, cfd.ied + 1):
            j = i - cfd.ist
            f.write(f"{cfd.mesh.xcc[j]:20.10e}{cfd.u[i]:20.10e}\n")
            
    # 可视化
    u_numerical = np.copy(cfd.u[cfd.ist:cfd.ied+1])
    print(f'u_numerical.size={u_numerical.size}')
    print(f'cfd.mesh.xcc.size={cfd.mesh.xcc.size}')
    #计算理论解
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
    
def visualizeEno(cfd):
    with open('solution.plt', 'w') as f:
        for i in range(cfd.ist, cfd.ied + 1):
            j = i - cfd.ist
            f.write(f"{cfd.mesh.xcc[j]:20.10e}{cfd.u[i]:20.10e}\n")
            
    # 可视化
    u_numerical = np.copy(cfd.u[cfd.ist:cfd.ied+1])
    print(f'u_numerical.size={u_numerical.size}')
    print(f'cfd.mesh.xcc.size={cfd.mesh.xcc.size}')
    #计算理论解
    u_analytical = analytical_solution(cfd.mesh.xcc, cfd.solver.T, cfd.solver.c, cfd.mesh.L)
    print(f'u_analytical.size={u_analytical.size}')
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.scatter(cfd.mesh.xcc, u_numerical, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=f'Numerical (Rusanov)')
    plt.plot(cfd.mesh.xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    ordinal1 = get_ordinal_numbers(cfd.iorder)
    ordinal2 = get_ordinal_numbers(cfd.solver.rk)
    plt.title(f'1D Convection Equation at t = {cfd.solver.T:.3f} using {cfd.iorder}{ordinal1}-order ENO and {cfd.solver.rk}{ordinal2}-order Runge-Kutta methods')
    plt.legend()
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)    
    plt.tight_layout()
    plt.show()
    
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
        
        for i in range(self.nnodes):
            self.x[i] = self.xmin + i * self.dx
            
        for i in range(self.ncells):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i+1])

class Solver:
    def __init__(self):
        #interpolation :0 Eno; 1 Weno
        self.interpolation = 0
        self.iflux = 0
        self.rk = 1
        self.c = 1.0
        self.T = 0.625
        #self.dt = .0025
        self.dt = .025
        
class Cfd:
    def __init__(self, solver, mesh, iorder):
        self.solver = solver
        self.mesh = mesh
        self.nx = mesh.nx
        self.iorder = iorder
        self.ighost = iorder
        self.ishift = self.ighost + 1
        self.ist = 0 + self.ishift
        self.ied = self.nx - 1 + self.ishift
        self.ntcell = self.nx + 2 * self.ishift
        self.isize = iorder * (iorder + 1)

        self.il = np.zeros(self.nx + 1, dtype=int)
        self.ir = np.zeros(self.nx + 1, dtype=int)
        self.coef = np.zeros((iorder + 1, iorder))
        self.dd = np.zeros((iorder, self.ntcell))
        self.up1_2m = np.zeros(self.nx + 1)
        self.up1_2p = np.zeros(self.nx + 1)
        self.flux = np.zeros(self.nx + 1)
        self.res = np.zeros(self.nx)

        # Field module variables
        self.u = np.zeros(self.ntcell)
        self.un = np.zeros(self.ntcell)
        
def RunEno(solver, mesh, iorder):
    cfd = Cfd(solver, mesh, iorder)
    init_coef(cfd.iorder, cfd.coef)
    init_field(cfd)
    
    simu_time = solver.T
    t = 0.0
    dt = solver.dt
    while t < simu_time:
        runge_kutta(cfd)
        if t + dt > simu_time:
            dt = simu_time - t
        t += dt
    print(f'T={t:.3f},runge_kutta{solver.rk},cfd{cfd.iorder},{fluxnames[solver.iflux]} FLUX')

    return np.copy(cfd.u[cfd.ist:cfd.ied+1])
    
def RunWeno(solver, mesh, iorder):
    cfd = Cfd(solver, mesh, iorder)
    init_coef(cfd.iorder, cfd.coef)
    init_field(cfd)
    
    simu_time = solver.T
    t = 0.0
    dt = solver.dt
    while t < simu_time:
        runge_kutta(cfd)
        if t + dt > simu_time:
            dt = simu_time - t
        t += dt
    print(f'T={t:.3f},runge_kutta{solver.rk},WENO{cfd.iorder},{fluxnames[solver.iflux]} FLUX')
    return np.copy(cfd.u[cfd.ist:cfd.ied+1])    
    
def performEnoOrderAnalysis():
    iorder_max = 7
    mesh = Mesh()
    solver = Solver()
    #计算理论解
    u_analytical = analytical_solution(mesh.xcc, solver.T, solver.c, mesh.L)
    
    u_list = []
    for iorder in range(1, iorder_max+1):
        u = RunEno(solver, mesh, iorder)
        u_list.append(u)
        
    plot_eno_OrderAnalysis(solver, mesh.xcc, u_list, u_analytical)
    
def performEnoTimestepAnalysis():
    mesh = Mesh()
    solver = Solver()
    #计算理论解
    u_analytical = analytical_solution(mesh.xcc, solver.T, solver.c, mesh.L)
    
    u_list = []
    dt_list = []
    solver.dt = 0.025/4
    #solver.dt = 0.025/16
    n = 12
    solver.rk = 3
    iorder = 7
    for i in range(0, n):
        u = RunEno(solver, mesh, iorder)
        u_list.append(u)
        dt_list.append(solver.dt)
        print(f'i={i+1},N={n},T={solver.T:.3f},dt={solver.dt},nt={int(solver.T/solver.dt)},runge_kutta{solver.rk},ENO{iorder},{fluxnames[solver.iflux]} FLUX')
        solver.dt /= 2
        
    plot_eno_TimestepAnalysis(solver, mesh.xcc, u_list, u_analytical, dt_list, iorder)
    
def performEnoWenoAnalysis():
    mesh = Mesh()
    solver = Solver()
    #计算理论解
    u_analytical = analytical_solution(mesh.xcc, solver.T, solver.c, mesh.L)
    
    solver.rk = 1
    solver.dt = 0.0025
    solver.iorder = 3
    
    u_list = []
    solver.interpolation = 0
    u = RunEno(solver, mesh, solver.iorder)
    u_list.append(u)
    solver.interpolation = 1
    u = RunWeno(solver, mesh, solver.iorder)
    u_list.append(u)
        
    plot_EnoWeno_Analysis(solver, mesh.xcc, u_list, u_analytical)
    
def plot_eno_OrderAnalysis(solver, xcc, u_list, u_analytical):
    # 定义一个包含不同颜色、线形和标记的列表
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
    plt.title(f'1D Convection Equation at t = {solver.T:.3f} using [1-7]th-order ENO and {solver.rk}{ordinal}-order Runge-Kutta methods')
    for i in range(0, n):
        lable = 'Numerical (Rusanov)ENO' + str(i+1)
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
    
def plot_EnoWeno_Analysis(solver, xcc, u_list, u_analytical):
    # 定义一个包含不同颜色、线形和标记的列表
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
    
def main():
    performEnoWenoAnalysis()

if __name__ == "__main__":
    main()