import numpy as np
import matplotlib.pyplot as plt

def conjugate_gradients(A, b):
    n = len(b)
    
    x = np.zeros( n )
    
    i = 0
    imax = 10000000000
    eps = 0.0001
    r = b - A @ x
    p = r
    #deltanew = r.T @ r
    deltanew = np.dot(r,r)
    delta0 = deltanew
    while i < imax and deltanew > eps**2 * delta0:
        #alpha = deltanew / (p.T @ (A @ p))
        alpha = deltanew / np.dot(p, (A @ p))
        x = x + alpha * p
        r = b - A @ x
        deltaold = deltanew
        #deltanew = r.T @ r
        deltanew = np.dot(r,r)
        beta = float( deltanew ) / float( deltaold )
        p = r + beta * p
        print("i=",i,"deltanew=",deltanew)
        i += 1
    print("niter=",i)
        
    return x    

n = 4
A = np.zeros( ( n, n ) )

for j in range(0, n):
    for i in range(0, n):
      if ( i== j ):
          if ( i!= (n-1) ):
              A[i,j] = 2
              A[i,j+1] = -1
              A[i+1,j] = -1
          else:
              A[i,j] = 1
              
L=1.0
nPoints = n+2
dx = L/(nPoints-1)
dx2= dx*dx
b = np.zeros( n )

for i in range(0, n-1):
    b[i] = dx2

b[n-1]=1.5*dx2

x = np.zeros( n )

for i in range(0, n):
    x[i] = (i+1)*dx

m = 100
f = np.zeros( m )
xx = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xx[i] - 0.5*xx[i]*xx[i]
    
   
#u = np.linalg.solve(A, b)
u = conjugate_gradients(A, b)

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
