import numpy as np
import matplotlib.pyplot as plt

def calc_ap( u, ap, dx, dy, nx, ny ):
    dx2 = dx * dx
    dy2 = dy * dy
    
    #boundary
    for j in range(0, ny):
        u[0,j] = 0.0
        u[nx-1,j] = u[nx-2,j]
        
    for i in range(0, nx):
        u[i,0] = u[i,1]
        u[i,ny-1] = u[i,ny-2]        
    
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            d2udx2 = ( u[i+1,j] - 2*u[i,j] + u[i-1,j] ) / dx2
            d2udy2 = ( u[i,j+1] - 2*u[i,j] + u[i,j-1] ) / dy2            
            ap[i,j] = d2udx2 + d2udy2
    return
    
def calc_residual( x, b, nx, ny, dx, dy ):
    ax = np.zeros( ( nx, ny ) )
    calc_ap( x, ax, dx, dy, nx, ny )
    r = b - ax
    return r
    
def my_dot_product( a, b, nx, ny ):
    sum = 0.0
    for j in range(1, ny-1):
        for i in range(1, nx-1):
            sum += a[i,j] * b[i,j]
    return sum

def conjugate_gradients( dx, dy,  b, nx, ny ):
    x = np.zeros( ( nx, ny ) )
    ap = np.zeros( ( nx, ny )  )
    r = np.zeros( ( nx, ny ) )
    
    i = 0
    imax = 10000000000
    eps = 0.0001
    
    r = calc_residual( x, b, nx, ny, dx, dy )
    p = r.copy()
    deltanew = my_dot_product( r, r, nx, ny )
    delta0 = deltanew

    while i < imax and deltanew > eps**2 * delta0:
        calc_ap( p, ap, dx, dy, nx, ny )
        alpha = deltanew / my_dot_product( p, ap, nx, ny )
        
        # update the numerical solution by adding some component of conjugate vector
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                x[i,j] = x[i,j] + alpha * p[i,j]
            
        # update the residual by removing some component of previous residual
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                r[i,j] = r[i,j] - alpha * ap[i,j]
       
        deltaold = deltanew
        deltanew = my_dot_product( r, r, nx, ny )
        beta = deltanew / deltaold
       
        # update the conjugate vector
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                p[i,j] = r[i,j] + beta * p[i,j]
        
        print("i=",i,"deltanew=",deltanew)
        i += 1
    print("niter=",i)
        
    return x    

n = 8
nx = n + 2
ny = n + 2
             
Lx = 1.0
dx = Lx/( nx - 1 )
dx2= dx * dx

Ly = 1.0
dy = Ly/( ny - 1 )
dy2= dy * dy

b = np.zeros( ( nx, ny ) )
for j in range(0, ny):
    for i in range(1, nx-2):
        b[i,j] = -1.0
    b[nx-2,j]=-1.5*1.0

#print("b=\n",b)

x = np.zeros( nx )
y = np.zeros( ny )

for i in range(0, nx):
    x[i] = i*dx

for j in range(0, ny):
    y[i] = j*dy
    
m = 100
f = np.zeros( m )
xx = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xx[i] - 0.5*xx[i]*xx[i]
    
   
#u = np.linalg.solve(A, b)
u = conjugate_gradients( dx, dy,  b, nx, ny )

if ( n <= 100 ):
    skip=1
else:
    skip=max(1,int(n/25))
    
mylabel="CG Calc N = {value}".format( value = n )
    
plt.figure("OneFLOW-CFD Solver Conjugate Gradient Methods", figsize=(6, 4), dpi=100)   
for k in range(0, ny):
    plt.plot(x[1:-1],u[1:-1,k],'bo', markersize=6,markerfacecolor='none',markevery=skip,label=mylabel )

plt.plot(xx,f,'k',linewidth=1,label="theory")
plt.legend()
plt.show()
