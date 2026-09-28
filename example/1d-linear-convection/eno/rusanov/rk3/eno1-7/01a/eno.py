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

def residual(q, eno):
    reconstruction(q, eno)
    inviscid_flux(eno.up1_2m, eno.up1_2p, eno.flux, eno)
    for i in range(eno.nx):
        eno.res[i] = -(eno.flux[i + 1] - eno.flux[i]) / eno.dx

def reconstruction(q, eno):
    # Choose the stencil by ENO method
    eno.dd[0, 0:eno.ntcell-1] = q[0:eno.ntcell-1]
    
    for m in range(1, eno.iorder):
        for j in range(0, eno.ntcell-1):
            eno.dd[m, j] = eno.dd[m-1, j+1] - eno.dd[m-1, j]
            
    for i in range(eno.nx + 1):
        eno.il[i] = i - 1
        for m in range(1, eno.iorder):
            if abs(dd[m, eno.il[i]-1+eno.ishift]) <= abs(dd[m, eno.il[i]+eno.ishift]):
                eno.il[i] -= 1
                
    for i in range(eno.nx + 1):
        eno.ir[i] = i
        for m in range(1, eno.iorder):
            if abs(dd[m, eno.ir[i]-1+eno.ishift]) <= abs(dd[m, ir[i]+eno.ishift]):
                eno.ir[i] -= 1
    
    # Reconstruction u(j+1/2)
    for i in range(eno.nx + 1):
        k1 = eno.il[i]
        k2 = eno.ir[i]
        l1 = i - k1
        l2 = i - k2
        eno.up1_2m[i] = 0
        eno.up1_2p[i] = 0
        for m in range(eno.iorder):
            eno.up1_2m[i] += q[k1 + eno.ishift + m] * eno.coef[l1, m]
            eno.up1_2p[i] += q[k2 + eno.ishift + m] * eno.coef[l2, m]
            
def inviscid_flux(up1_2m, up1_2p, flux, eno):
    iflux = 0
    if iflux == 0:
        rusanov_flux(up1_2m, up1_2p, flux, eno)
    else:
        engquist_osher_flux(up1_2m, up1_2p, flux, eno)

def engquist_osher_flux(up1_2m, up1_2p, flux, eno):
    for i in range(eno.nx + 1):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        
        cp = 0.5 * ( eno.c + abs(eno.c) )
        cm = 0.5 * ( eno.c - abs(eno.c) )
        
        flux[i] = cp * uL + cm * u_R
        
def rusanov_flux(up1_2m, up1_2p, flux, eno):
    for i in range(eno.nx + 1):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        F_L = eno.c * u_L  # 左状态通量
        F_R = eno.c * u_R  # 右状态通量
        alpha = abs(eno.c) # 最大波速
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)
                

def boundary(u, eno):
    for i in range(-eno.ighost, 1):
        u[eno.ist - 1 + i] = u[eno.ied + i]
    for i in range(1, eno.ighost + 2):
        u[eno.ied + i] = u[eno.ist - 1 + i]

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

def init_field(eno):
    for i in range(eno.ist, eno.ied + 1):
        j = i - eno.ist
        if 0.5 <= eno.xcc[j] <= 1.0:
            eno.u[i] = 2.0
        else:
            eno.u[i] = 1.0
    boundary(eno.u, eno)
    update_oldfield(eno.un, eno.u)
    
def runge_kutta(eno):
    rk = eno.rk
    if rk == 1:
        runge_kutta_1(eno)
    elif rk == 2:
        runge_kutta_2(eno)
    else:
        runge_kutta_3(eno)
    
def runge_kutta_1(eno):
    global u, un, dt
    residual(u, eno)
    for i in range(nx):
        j = i + ishift
        u[j] = u[j] + dt * res[i]
    boundary(u)
    update_oldfield(un, u)
    
def runge_kutta_2(eno):
    #global u, un, dt
    residual(eno.u, eno)
    for i in range(eno.nx):
        j = i + eno.ishift
        eno.u[j] = eno.u[j] + eno.dt * eno.res[i]
    boundary(eno.u, eno)
    
    residual(eno.u, eno)
    for i in range(eno.nx):
        j = i + eno.ishift
        eno.u[j] = 0.5 * eno.un[j] + 0.5 * eno.u[j] + 0.5 * eno.dt * eno.res[i]
    boundary(eno.u, eno)
    update_oldfield(eno.un, eno.u)

def runge_kutta_3(eno):
    global u, un, dt
    residual(eno.u, eno)
    for i in range(eno.nx):
        j = i + eno.ishift
        eno.u[j] = eno.u[j] + eno.dt * eno.res[i]
    boundary(eno.u, eno)
    
    residual(eno.u, eno)
    for i in range(eno.nx):
        j = i + eno.ishift
        eno.u[j] = 0.75 * eno.un[j] + 0.25 * eno.u[j] + 0.25 * eno.dt * eno.res[i]
    boundary(eno.u, eno)
    
    residual(eno.u, eno)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in range(eno.nx):
        j = i + eno.ishift
        eno.u[j] = c1 * eno.un[j] + c2 * eno.u[j] + c3 * eno.dt * eno.res[i]
    boundary(eno.u, eno)
    update_oldfield(eno.un, eno.u)

def get_ordinal_numbers(order):
    if order == 1:
        return 'st'
    elif order == 2:
        return 'nd'
    elif order == 3:
        return 'rd'
    else:
        return 'th'

def visualize(eno):
    with open('solution.plt', 'w') as f:
        for i in range(eno.ist, eno.ied + 1):
            j = i - eno.ist
            f.write(f"{eno.xcc[j]:20.10e}{eno.u[i]:20.10e}\n")
            
    # 可视化
    u_numerical = np.copy(eno.u[eno.ist:eno.ied+1])
    print(f'u_numerical.size={u_numerical.size}')
    print(f'eno.xcc.size={eno.xcc.size}')
    #计算理论解
    u_analytical = analytical_solution(eno.xcc, eno.T, eno.c, eno.L)    
    print(f'u_analytical.size={u_analytical.size}')
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.scatter(eno.xcc, u_numerical, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=f'Numerical (Rusanov)')
    plt.plot(eno.xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    ordinal1 = get_ordinal_numbers(eno.iorder)
    ordinal2 = get_ordinal_numbers(eno.rk)
    plt.title(f'1D Convection Equation at t = {eno.T:.3f} using {eno.iorder}{ordinal1}-order ENO and {eno.rk}{ordinal2}-order Runge-Kutta methods')
    plt.legend()
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)    
    plt.tight_layout()
    plt.show()
            
class Eno:
    def __init__(self, iorder):
        self.nx = 40
        self.rk = 2
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
        self.L = 2.0
        self.x = np.zeros(self.nx + 1)
        self.xcc = np.zeros(self.nx)

        # Field module variables
        self.u = np.zeros(self.ntcell)
        self.un = np.zeros(self.ntcell)
        
        self.c = 1.0
        self.T = 0.625
        self.dt = .0025
    def init_mesh(self):
        xstart = 0.0
        xend   = self.L
        self.dx = (xend - xstart) / self.nx
        
        for i in range(0, self.nx+1):
            self.x[i] = xstart + i * self.dx
            
        for i in range(0, self.nx):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i + 1])
            
def RunEno(iorder):
    eno = Eno(iorder)
    init_coef(eno.iorder, eno.coef)
    eno.init_mesh()
    init_field(eno)
    
    nt = int(eno.T / eno.dt)
    simu_time = eno.T
    
    t = 0.0
    
    while t < simu_time:
        runge_kutta(eno)
        t += eno.dt
        if t + eno.dt > simu_time:
            dt = simu_time - t
    
    print(t)
    visualize(eno)
    
def main():
    iorder = 1
    RunEno(iorder)

if __name__ == "__main__":
    main()