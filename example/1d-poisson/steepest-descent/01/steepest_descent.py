import numpy as np
import matplotlib.pyplot as plt

def steepest_descent(A, b):
    n = len(b)
    print("len(b)=",n)
    print("b=\n",b)
    
    x = np.zeros( n )
    
    i = 0
    imax = 1000
    eps = 0.01
    #print("A * x=\n",A * x)
    #print("A @ x=\n",A @ x)
    #r = b - A * x
    r = b - A @ x
    #print("r=\n",r)
    #print("x=\n",x)
    delta = r.T @ r
    delta0 = delta
    #print("delta0=\n",delta0)
    #print("delta=\n",delta)
        
    while i < imax and delta > eps**2 * delta0:
        term2 = r.T @ (A @ r)
        alpha = delta / term2
        #print("delta=\n",delta)
        #print("term2=\n",term2)
        #print("alpha=\n",alpha)
        x = x + alpha * r
        #print("x=\n",x)
        r = b - A @ x
        delta = r.T @ r
        print("i=",i,"delta=",delta)
        i += 1    
    print("niter=",i)
    return n

n = 400
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
              
#print("A=\n",A)

L=1.0
nPoints = n+2
dx = L/(nPoints-1)
dx2= dx*dx
b = np.zeros( n )

for i in range(0, n-1):
    b[i] = dx2

b[n-1]=0.5*dx2

x = np.zeros( n )

for i in range(0, n):
    x[i] = (i+1)*dx

m = 10
f = np.zeros( m )
xx = np.linspace(0,1,m)

for i in range(0, m):
    f[i] = xx[i] - 0.5*xx[i]*xx[i]
    
u = np.linalg.solve(A, b)

#u  = steepest_descent(A, b)
steepest_descent(A, b)

plt.plot(x,u,'bo', markersize=6,markevery=10,label="OneFLOW-CFD" )
plt.plot(xx,f,'k',linewidth=1,label="theory")
plt.legend()
plt.show()
