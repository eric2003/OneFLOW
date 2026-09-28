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
def analytical_solution(x, t, a):
    # 初始条件沿 x - at 平移
    x_shifted = x - a * t
    return initial_condition((x_shifted + L) % L)  # 周期边界条件  

def residual(q):
    reconstruction(q)
    inviscid_flux(up1_2m, up1_2p, flux)
    for i in range(nx):
        res[i] = -(flux[i + 1] - flux[i]) / dx

def reconstruction(q):
    global il, ir, dd, up1_2m, up1_2p
    
    # Choose the stencil by ENO method
    dd[0, 0:ntcell-1] = q[0:ntcell-1]
    
    for m in range(1, iorder):
        for j in range(0, ntcell-1):
            dd[m, j] = dd[m-1, j+1] - dd[m-1, j]
            
    for i in range(nx + 1):
        il[i] = i - 1
        for m in range(1, iorder):
            if abs(dd[m, il[i]-1+ishift]) <= abs(dd[m, il[i]+ishift]):
                il[i] -= 1
                
    for i in range(nx + 1):
        ir[i] = i
        for m in range(1, iorder):
            if abs(dd[m, ir[i]-1+ishift]) <= abs(dd[m, ir[i]+ishift]):
                ir[i] -= 1
    
    # Reconstruction u(j+1/2)
    for i in range(nx + 1):
        k1 = il[i]
        k2 = ir[i]
        l1 = i - k1
        l2 = i - k2
        up1_2m[i] = 0
        up1_2p[i] = 0
        for m in range(iorder):
            up1_2m[i] += q[k1 + ishift + m] * coef[l1, m]
            up1_2p[i] += q[k2 + ishift + m] * coef[l2, m]
            
def inviscid_flux(up1_2m, up1_2p, flux):
    iflux = 0
    if iflux == 0:
        rusanov_flux(up1_2m, up1_2p, flux)
    else:
        engquist_osher_flux(up1_2m, up1_2p, flux)

def engquist_osher_flux(up1_2m, up1_2p, flux):
    for i in range(nx + 1):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        
        cp = 0.5 * ( c + abs(c) )
        cm = 0.5 * ( c - abs(c) )
        
        flux[i] = cp * uL + cm * u_R
        
def rusanov_flux(up1_2m, up1_2p, flux):
    for i in range(nx + 1):
        u_L = up1_2m[i]
        u_R = up1_2p[i]
        F_L = c * u_L  # 左状态通量
        F_R = c * u_R  # 右状态通量
        alpha = abs(c) # 最大波速
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)
                

def boundary(u):
    for i in range(-ighost, 1):
        u[ist - 1 + i] = u[ied + i]
    for i in range(1, ighost + 2):
        u[ied + i] = u[ist - 1 + i]

def update_oldfield(qn, q):
    qn[:] = q[:]

def init_coef():
    global coef
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

def init_mesh():
    global xstart, xend, dx, x, xcc
    xstart = 0.0
    xend   = L
    dx = (xend - xstart) / nx
    
    for i in range(0, nx+1):
        x[i] = xstart + i * dx
        
    for i in range(0, nx):
        xcc[i] = 0.5 * (x[i] + x[i + 1])

def init_field():
    global u, un
    for i in range(ist, ied + 1):
        j = i - ist
        if 0.5 <= xcc[j] <= 1.0:
            u[i] = 2.0
        else:
            u[i] = 1.0
    boundary(u)
    update_oldfield(un, u)
    
def runge_kutta():
    if rk == 1:
        runge_kutta_1()
    elif rk == 2:
        runge_kutta_2()
    else:
        runge_kutta_3()
    
def runge_kutta_1():
    global u, un, dt
    residual(u)
    for i in range(nx):
        j = i + ishift
        u[j] = u[j] + dt * res[i]
    boundary(u)
    update_oldfield(un, u)
    
def runge_kutta_2():
    global u, un, dt
    residual(u)
    for i in range(nx):
        j = i + ishift
        u[j] = u[j] + dt * res[i]
    boundary(u)
    
    residual(u)
    for i in range(nx):
        j = i + ishift
        u[j] = 0.5 * un[j] + 0.5 * u[j] + 0.5 * dt * res[i]
    boundary(u)    
    update_oldfield(un, u)

def runge_kutta_3():
    global u, un, dt
    residual(u)
    for i in range(nx):
        j = i + ishift
        u[j] = u[j] + dt * res[i]
    boundary(u)
    
    residual(u)
    for i in range(nx):
        j = i + ishift
        u[j] = 0.75 * un[j] + 0.25 * u[j] + 0.25 * dt * res[i]
    boundary(u)
    
    residual(u)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in range(nx):
        j = i + ishift
        u[j] = c1 * un[j] + c2 * u[j] + c3 * dt * res[i]
    boundary(u)
    update_oldfield(un, u)

def get_ordinal_numbers(order):
    if order == 1:
        return 'st'
    elif order == 2:
        return 'nd'
    elif order == 3:
        return 'rd'
    else:
        return 'th'

def visualize():
    with open('solution.plt', 'w') as f:
        for i in range(ist, ied + 1):
            j = i - ist
            f.write(f"{xcc[j]:20.10e}{u[i]:20.10e}\n")
            
    # 可视化
    u_numerical = np.copy(u[ist:ied+1])
    print(f'u_numerical.size={u_numerical.size}')
    print(f'xcc.size={xcc.size}')
    #计算理论解
    u_analytical = analytical_solution(xcc, T, c)    
    print(f'u_analytical.size={u_analytical.size}')
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.scatter(xcc, u_numerical, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=f'Numerical (Rusanov)')
    plt.plot(xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    ordinal1 = get_ordinal_numbers(iorder)
    ordinal2 = get_ordinal_numbers(rk)
    plt.title(f'1D Convection Equation at t = {T:.3f} using {iorder}{ordinal1}-order ENO and {rk}{ordinal2}-order Runge-Kutta methods')
    plt.legend()
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)    
    plt.tight_layout()
    plt.show()            
            
# Global constants and variables
nx = 40
rk = 2
iorder = 7
ighost = iorder
ishift = ighost + 1
ist = 0 + ishift
ied = nx - 1 + ishift
ntcell = nx + 2 * ishift
isize = iorder * (iorder + 1)

il = np.zeros(nx + 1, dtype=int)
ir = np.zeros(nx + 1, dtype=int)
coef = np.zeros((iorder + 1, iorder))
dd = np.zeros((iorder, ntcell))
up1_2m = np.zeros(nx + 1)
up1_2p = np.zeros(nx + 1)
flux = np.zeros(nx + 1)
res = np.zeros(nx)
dt = 0.0

# Mesh module variables
xstart = 0.0
xend = 0.0
L = 2.0
dx = 0.0
x = np.zeros(nx + 1)
xcc = np.zeros(nx)

# Field module variables
u = np.zeros(ntcell)
un = np.zeros(ntcell)

def main():
    global dt, T, c
    init_coef()
    init_mesh()
    init_field()
    
    c = 1.0
    T = 0.625
    dt = .0025
    nt = int(T / dt)
    simu_time = T
    
    t = 0.0
    
    while t < simu_time:
        runge_kutta()
        t += dt
        if t + dt > simu_time:
            dt = simu_time - t
    
    print(t)
    visualize()

if __name__ == "__main__":
    main()