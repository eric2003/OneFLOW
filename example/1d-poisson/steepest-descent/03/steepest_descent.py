import numpy as np
import matplotlib.pyplot as plt

def steepest_descent(A, b):
    n = len(b)
    #print("len(b)=",n)
    #print("b=\n",b)
    
    x = np.zeros( n )
    
    i = 0
    imax = 10000000000
    eps = 0.01
    r = b - A @ x
    delta = r.T @ r
    delta0 = delta
    #epsl=1.0e-12
    #epsl=1.0e-13
    #epsl=1.0e-14
    epsl=1.0e-15
        
    #while i < imax and delta > eps**2 * delta0:
    while i < imax and delta > epsl:
        alpha = delta / ( r.T @ (A @ r) )
        x = x + alpha * r
        r = b - A @ x
        delta = r.T @ r
        print("i=",i,"delta=",delta)
        i += 1    
    print("niter=",i)
    return x

n = 128

L=1.0
nPoints = n+2
dx = L/(nPoints-1)
dx2= dx*dx

A = np.zeros( ( n, n ) )
for j in range(0, n):
    for i in range(0, n):
      if ( i== j ):
          if ( i!= (n-1) ):
              A[i,j] = 2.0/dx2
              A[i,j+1] = -1.0/dx2
              A[i+1,j] = -1.0/dx2
          else:
              A[i,j] = 1.0/dx2
              

b = np.zeros( n )

for i in range(0, n-1):
    b[i] = 1.0

b[n-1]=1.5*1.0

x = np.zeros( n )

for i in range(0, n):
    x[i] = (i+1)*dx

m = 100
f = np.zeros( m )
xx = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xx[i] - 0.5*xx[i]*xx[i]
    
   
#u = np.linalg.solve(A, b)
u = steepest_descent(A, b)
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
