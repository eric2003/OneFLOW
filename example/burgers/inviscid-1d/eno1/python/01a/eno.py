import numpy as np

def residual(q):
    reconstruction(q)
    engquist_osher_flux(up1_2m, up1_2p, flux)
    for i in range(nx):
        res[i] = -(flux[i + 1] - flux[i]) / dx

def reconstruction(q):
    global il, ir, dd, up1_2m, up1_2p
    
    # Choose the stencil by ENO method
    dd[0, 1:ntcell + 1] = q[1:ntcell + 1]
    
    for m in range(1, iorder):
        for j in range(1, ntcell):
            dd[m, j] = dd[m-1, j+1] - dd[m-1, j]
    
    for i in range(nx + 1):
        il[i] = i
        ir[i] = i + 1
        for m in range(1, iorder):
            if abs(dd[m, il[i]-1+ishift]) <= abs(dd[m, il[i]+ishift]):
                il[i] -= 1
            if abs(dd[m, ir[i]-1+ishift]) <= abs(dd[m, ir[i]+ishift]):
                ir[i] -= 1
    
    # Reconstruction u(j+1/2)
    for i in range(nx + 1):
        k1 = il[i]
        k2 = ir[i]
        l1 = i - k1 + 1
        l2 = i - k2 + 1
        up1_2m[i] = 0
        up1_2p[i] = 0
        for m in range(iorder):
            up1_2m[i] += q[k1 + ishift + m] * coef[l1, m]
            up1_2p[i] += q[k2 + ishift + m] * coef[l2, m]

def engquist_osher_flux(up1_2m, up1_2p, flux):
    for i in range(nx + 1):
        if up1_2m[i] >= 0:
            if up1_2p[i] >= 0:
                flux[i] = 0.5 * up1_2m[i] * up1_2m[i]
            else:
                flux[i] = 0.5 * (up1_2m[i] * up1_2m[i] + up1_2p[i] * up1_2p[i])
        else:
            if up1_2p[i] >= 0:
                flux[i] = 0
            else:
                flux[i] = 0.5 * up1_2p[i] * up1_2p[i]

def boundary(u):
    for i in range(-ighost, 1):
        u[ishift + i] = u[ied + i]
    for i in range(1, ighost + 1):
        u[ied + i] = u[ishift + i]

def update_oldfield(qn, q):
    qn[:] = q[:]

def init_coef():
    global coef
    coef[0] = [1.0]
    coef[1] = [1.0]
def init_mesh():
    global xstart, xend, dx, x, xcc
    xstart = -1.0
    xend = 1.0
    dx = (xend - xstart) / nx
    xstart0 = xstart - ishift * dx
    
    for i in range(1, ntcell + 2):
        x[i] = xstart0 + (i - 1) * dx
    
    for i in range(1, ntcell + 1):
        xcc[i] = 0.5 * (x[i] + x[i + 1])

def init_field():
    global u, un
    for i in range(ist, ied + 1):
        u[i] = 0.25 + 0.5 * np.sin(pi * xcc[i])
    boundary(u)
    update_oldfield(un, u)

def runge_kutta_3():
    global u, un, dt
    residual(u)
    for i in range(nx):
        j = i + 1 + ishift
        u[j] = u[j] + dt * res[i]
    boundary(u)
    
    residual(u)
    for i in range(nx):
        j = i + 1 + ishift
        u[j] = 0.75 * un[j] + 0.25 * u[j] + 0.25 * dt * res[i]
    boundary(u)
    
    residual(u)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in range(nx):
        j = i + 1 + ishift
        u[j] = c1 * un[j] + c2 * u[j] + c3 * dt * res[i]
    boundary(u)
    update_oldfield(un, u)

def visualize():
    with open('solution_total.plt', 'w') as f1:
        for i in range(1, ntcell + 1):
            f1.write(f"{xcc[i]:20.10e}{u[i]:20.10e}\n")
    
    with open('solution.plt', 'w') as f2:
        for i in range(ist, ied + 1):
            f2.write(f"{xcc[i]:20.10e}{u[i]:20.10e}\n")
            
# Global constants and variables
nx = 40
iorder = 1
ighost = iorder
ishift = ighost + 1
ist = 1 + ishift
ied = nx + ishift
ntcell = nx + ishift + ighost
isize = iorder * (iorder + 1)
pi = 3.14159265358979323846

il = np.zeros(nx + 1, dtype=int)
ir = np.zeros(nx + 1, dtype=int)
coef = np.zeros((iorder + 1, iorder))
dd = np.zeros((ighost, ntcell + 1))
up1_2m = np.zeros(nx + 1)
up1_2p = np.zeros(nx + 1)
flux = np.zeros(nx + 1)
res = np.zeros(nx)
dt = 0.0

# Mesh module variables
xstart = 0.0
xend = 0.0
dx = 0.0
x = np.zeros(ntcell + 2)
xcc = np.zeros(ntcell + 1)

# Field module variables
u = np.zeros(ntcell + 1)
un = np.zeros(ntcell + 1)

def main():
    global dt
    init_coef()
    init_mesh()
    init_field()
    
    simu_time = float(input("Input T: "))
    dt = dx * 0.5
    t = 0.0
    
    while t < simu_time:
        runge_kutta_3()
        t += dt
        if t + dt > simu_time:
            dt = simu_time - t
    
    print(t)
    visualize()

if __name__ == "__main__":
    main()