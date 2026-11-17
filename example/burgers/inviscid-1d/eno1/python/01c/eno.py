import numpy as np

def residual(q):
    reconstruction(q)
    engquist_osher_flux(up1_2m, up1_2p, flux)
    for i in range(nx):
        res[i] = -(flux[i + 1] - flux[i]) / dx

def reconstruction(q):
    global il, ir, dd, up1_2m, up1_2p
    
    # Choose the stencil by ENO method
    #dd[0, 1:ntcell + 1] = q[1:ntcell + 1]
    dd[0, 0:ntcell-1] = q[0:ntcell-1]
    
    for m in range(1, iorder):
        #for j in range(1, ntcell):
        for j in range(0, ntcell-1):
            dd[m, j] = dd[m-1, j+1] - dd[m-1, j]
            
    #ighost = 2
    #ishift = ighost + 1=3
    #i=0,1,...,nx
    #i+ishift=3,4,...,nx+3
   
    for i in range(nx + 1):
        il[i] = i - 1
        for m in range(1, iorder):
            if abs(dd[m, il[i]-1+ishift]) <= abs(dd[m, il[i]+ishift]):
                il[i] -= 1
                
    #i=0,1,...,nx
    #i+ishift=3,4,...,nx+3
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
        #print(f'l1={l1},l2={l2},k1={k1},k2={k2},i={i}')
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
    #ighost=2,-ighost=-2
    #ishift = ighost + 1=3
    #ist = 1 + ishift
    #ied = nx + ishift    
    #i=-2,-1,0
    #ishift+i=3-2,3-1,3+0=1,2,3
    #ied + i=ied-2,ied-1,ied+0
    #ist=ishift
    #ied=nx-1+ishift
    #ist-1=2
    for i in range(-ighost, 1):
        u[ist - 1 + i] = u[ied + i]
    #i=1,2
    #ied+i=ied+1,ied+2
    #ist-1=2
    #ist-1+i=3,4
    for i in range(1, ighost + 1):
        u[ied + i] = u[ist - 1 + i]

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
    
    for i in range(0, nx+1):
        x[i] = xstart + i * dx
        
    for i in range(0, nx):
        xcc[i] = 0.5 * (x[i] + x[i + 1])

def init_field():
    global u, un
    for i in range(ist, ied + 1):
        j = i - ist
        u[i] = 0.25 + 0.5 * np.sin(pi * xcc[j])
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

def visualize():
    with open('solution.plt', 'w') as f:
        for i in range(ist, ied + 1):
            j = i - ist
            f.write(f"{xcc[j]:20.10e}{u[i]:20.10e}\n")
            
# Global constants and variables
nx = 40
iorder = 1
#ighost = iorder
ighost = 2
ishift = ighost + 1
#ist = 1 + ishift
#ied = nx + ishift
ist = 0 + ishift
ied = nx - 1 + ishift
ntcell = nx + ishift + ighost
isize = iorder * (iorder + 1)
pi = 3.14159265358979323846

il = np.zeros(nx + 1, dtype=int)
ir = np.zeros(nx + 1, dtype=int)
coef = np.zeros((iorder + 1, iorder))
#dd = np.zeros((ighost, ntcell + 1))
dd = np.zeros((ighost, ntcell))
up1_2m = np.zeros(nx + 1)
up1_2p = np.zeros(nx + 1)
flux = np.zeros(nx + 1)
res = np.zeros(nx)
dt = 0.0

# Mesh module variables
xstart = 0.0
xend = 0.0
dx = 0.0
x = np.zeros(nx + 1)
xcc = np.zeros(nx)

# Field module variables
u = np.zeros(ntcell)
un = np.zeros(ntcell)

def main():
    global dt
    init_coef()
    init_mesh()
    init_field()
    
    simu_time = float(input("Input T: "))
    dt = dx * 0.5
    print(f'dt={dt}')
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