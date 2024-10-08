import numpy as np
import matplotlib.pyplot as plt

def calc_ap( u, ap, dx, n ):
    dx2 = dx * dx
    d2udx2 = ( - 2*u[0] + u[1] ) / dx2
    ap[0] = d2udx2
    #i=1,2,...,n-3,n-2
    for i in range(1, n-1):
        d2udx2 = ( u[i+1] - 2*u[i] + u[i-1] ) / dx2
        ap[i] = d2udx2
        
    d2udx2 = ( u[n-2] - u[n-1] ) / dx2
    ap[n-1] = d2udx2
    return            
    
def calc_residual( x, b, r, n, dx):
    ax = np.zeros( n )
    calc_ap( x, ax, dx, n )
    r[:] = b - ax    
    
def steepest_descent(dx, b):
    n = len(b)
    
    x = np.zeros( n )
    ax = np.zeros( n )
    r = np.zeros( n )
    
    i = 0
    imax = 10000000000
    eps = 0.01
    calc_residual( x, b, r, n, dx )
    print("r=\n",r)
    delta = np.dot(r, r)
    print("delta=\n",delta)
    delta0 = delta
    #epsl=1.0e-12
    #epsl=1.0e-13
    #epsl=1.0e-14
    epsl=1.0e-15
        
    while i < imax and delta > epsl:
        calc_ap( r, ax, dx, n )
        alpha = delta / np.dot(r, ax)
        x = x + alpha * r
        calc_residual( x, b, r, n, dx )
        delta = np.dot(r, r)
        print("i=",i,"delta=",delta)
        i += 1    
    print("niter=",i)
    return x

n = 4

L=1.0
nPoints = n+2
dx = L/(nPoints-1)
dx2= dx*dx


b = np.zeros( n )

for i in range(0, n-1):
    b[i] = -1.0

b[n-1]=-1.5*1.0

x = np.zeros( n )

for i in range(0, n):
    x[i] = (i+1)*dx

m = 100
f = np.zeros( m )
xx = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xx[i] - 0.5*xx[i]*xx[i]
    
   
#u = np.linalg.solve(A, b)
u = steepest_descent(dx, b)
if ( n <= 100 ):
    skip=1
else:
    skip=max(1,int(n/25))
    
mylabel="Calc N = {value}".format( value = n )
    
plt.figure("OneFLOW-CFD Solver Steepest Descent Methods", figsize=(6, 4), dpi=100)    
plt.plot(x,u,'bo', markersize=6,markerfacecolor='none',markevery=skip,label=mylabel )
plt.plot(xx,f,'k',linewidth=1,label="theory")
plt.legend()
plt.show()
