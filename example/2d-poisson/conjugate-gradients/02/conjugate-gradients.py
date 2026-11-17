import numpy as np
import matplotlib.pyplot as plt

def boundary( u, dx, dy, nx, ny ):
    dx2 = dx * dx
    dy2 = dy * dy
    
    #boundary
    for j in range(0, ny):
        u[0,j] = 0.0
        u[nx-1,j] = u[nx-2,j] + 0.5 * dx2
        
    for i in range(0, nx):
        u[i,0] = u[i,1]
        u[i,ny-1] = u[i,ny-2]
    return

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
    y[j] = j*dy
    
xx, yy = np.meshgrid(x, y, indexing='ij')

ue = np.zeros( ( nx, ny ) )
for j in range(0, ny):
    for i in range(0, nx):
        ue[i,j] = x[i] - 0.5 * x[i] * x[i]
    
m = 100
f = np.zeros( m )
xt = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xt[i] - 0.5*xt[i]*xt[i]
    
   
#u = np.linalg.solve(A, b)
u = conjugate_gradients( dx, dy, b, nx, ny )
boundary( u, dx, dy, nx, ny )

uerror = np.zeros( ( nx, ny ) )
uerror = u - ue

max_error = np.max( np.abs(uerror) )

#print("uerror=\n",uerror)        

print("Maximum Norm = ", max_error);


if ( n <= 100 ):
    skip=1
else:
    skip=max(1,int(n/25))
    
mylabel="CG Calc N = {value}".format( value = n )
    
fig = plt.figure("OneFLOW-CFD Solver+2D Poisson Equation+Conjugate Gradients", figsize=(10, 4), dpi=100)

ax1 = fig.add_subplot(1,3,1)
cs1 = ax1.contourf( xx, yy, ue, levels=20, cmap="jet",vmin=-1,vmax=1)
fig.colorbar(cs1, ax = ax1)
ax1.set_title("Exact solution")

ax2 = fig.add_subplot(1,3,2)
cs2 = ax2.contourf( xx, yy, u, levels=20, cmap="jet",vmin=-1,vmax=1)
fig.colorbar(cs2, ax = ax2)
ax2.set_title("Numerical solution")

fig.tight_layout()
fig.savefig("cg_contour.pdf")

ax3 = fig.add_subplot(1,3,3)

#for k in range(0, ny):
k=int(ny/2)
ax3.plot(x[1:-1],u[1:-1,k],'bo', markersize=6,markerfacecolor='none',markevery=skip,label=mylabel )

ax3.plot(xt,f,'k',linewidth=1,label="theory")
plt.legend()

plt.show()
