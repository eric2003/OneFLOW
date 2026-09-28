import numpy as np
import matplotlib.pyplot as plt

def calc_ap( u, ap, dx, nx ):
    dx2 = dx * dx
    d2udx2 = ( - 2*u[1] + u[2] ) / dx2
    ap[1] = d2udx2
    #i=1,2,...,nx-4,nx-3
    for i in range(1, nx-2):
        d2udx2 = ( u[i+1] - 2*u[i] + u[i-1] ) / dx2
        ap[i] = d2udx2
        
    d2udx2 = ( u[nx-3] - u[nx-2] ) / dx2
    ap[nx-2] = d2udx2
    return            
    
def calc_residual( x, b, nx, dx ):
    ax = np.zeros( nx )
    calc_ap( x, ax, dx, n )
    r = b - ax
    return r
    
def my_dot_product( a, b, nx ):
    sum = 0.0
    for i in range(1, nx-1):
        sum += a[i] * b[i]
    return sum    

def conjugate_gradients(dx, b, nx):
    x = np.zeros( nx  )
    ap = np.zeros( nx  )
    r = np.zeros( nx )
    
    i = 0
    imax = 10000000000
    eps = 0.0001
    
    r = calc_residual( x, b, nx, dx )
    p = r.copy()
    deltanew = my_dot_product(r,r,nx)
    delta0 = deltanew

    while i < imax and deltanew > eps**2 * delta0:
        calc_ap( p, ap, dx, nx )
        alpha = deltanew / my_dot_product( p, ap, nx )
        
        # update the numerical solution by adding some component of conjugate vector
        for i in range(1, nx-1):
            x[i] = x[i] + alpha * p[i]
            
        # update the residual by removing some component of previous residual
        for i in range(1, nx-1):
            r[i] = r[i] - alpha * ap[i]
       
        deltaold = deltanew
        deltanew = my_dot_product( r, r, nx )
        beta = deltanew / deltaold
       
        # update the conjugate vector
        for i in range(1, nx-1):
            p[i] = r[i] + beta * p[i]
        
        print("i=",i,"deltanew=",deltanew)
        i += 1
    print("niter=",i)
        
    return x    

n = 32
nx = n + 2
             
L=1.0
nPoints = n+2
dx = L/(nPoints-1)
dx2= dx*dx
b = np.zeros( nx )

for i in range(1, nx-2):
    b[i] = -1.0

b[nx-2]=-1.5*1.0

#print("b=\n",b)

x = np.zeros( nx )

for i in range(0, nx):
    x[i] = i*dx
    
#print("x=\n",x)
#exit()

m = 100
f = np.zeros( m )
xx = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xx[i] - 0.5*xx[i]*xx[i]
    
   
#u = np.linalg.solve(A, b)
u = conjugate_gradients(dx, b, nx)

if ( n <= 100 ):
    skip=1
else:
    skip=max(1,int(n/25))
    
mylabel="CG Calc N = {value}".format( value = n )
    
plt.figure("OneFLOW-CFD Solver Conjugate Gradients Methods", figsize=(6, 4), dpi=100)    
plt.plot(x,u,'bo', markersize=6,markerfacecolor='none',markevery=skip,label=mylabel )
plt.plot(xx,f,'k',linewidth=1,label="theory")
plt.legend()
plt.show()
