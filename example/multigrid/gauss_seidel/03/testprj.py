import numpy as np
import matplotlib.pyplot as plt

def Jacobi(v, niter, dx ):
    v2 = v.copy()
    coef = 1.0/3.0
    res = np.zeros(niter)
    for iter in range(0, niter):
        for i in range(1, n):
            v2[i] = coef * ( v[i-1] + v[i] + v[i+1] )
        v = v2.copy()
        resm = np.abs(np.max(v2))
        res[iter] = resm
    return res

def Gauss_Seidel(v, niter, dx):
    dx2 = dx * dx
    den = - 2.0 / dx2
    res = np.zeros(niter)
    for iter in range(0, niter):
        for i in range(1, n):
            d2udx2 = ( v[i-1] - 2*v[i] + v[i+1] )/dx2
            rm = -d2udx2
            v[i] = v[i] + rm / den
        resm = np.abs(np.max(v))
        res[iter] = resm
    return res
    
    
def Red_Black_Gauss_Seidel(v, niter, dx):
    dx2 = dx * dx
    den = - 2.0 / dx2
    res = np.zeros(niter)
    for iter in range(0, niter):
        #red sweep
        for i in range(1, n, 2):
            d2udx2 = ( v[i-1] - 2*v[i] + v[i+1] )/dx2
            rm = -d2udx2
            v[i] = v[i] + rm / den
        #black sweep
        for i in range(2, n, 2):
            d2udx2 = ( v[i-1] - 2*v[i] + v[i+1] )/dx2
            rm = -d2udx2
            v[i] = v[i] + rm / den
        resm = np.abs(np.max(v))
        res[iter] = resm
    return res    

n = 64
h = 1.0/n
dx = h

print("dx=",dx)

x = np.zeros( n+1 )
y1 = np.zeros( n+1 )
y3 = np.zeros( n+1 )
y6 = np.zeros( n+1 )
v = np.zeros( n+1 )

for i in range(0, n+1):
    x[i] = i * dx
    y1[i] = np.sin( i * 1 * np.pi / n )
    y3[i] = np.sin( i * 3 * np.pi / n )
    y6[i] = np.sin( i * 6 * np.pi / n )
    
#print("x=\n",x)


#print("v=\n",v)
#resv = np.abs(np.max(v))
#print("resv=",resv)

niter = 100

v = y1.copy()
res1Jacobi = Jacobi(v, niter, dx)

v = y1.copy()
res1GS = Gauss_Seidel(v, niter, dx)

v = y1.copy()
res1RBGS = Red_Black_Gauss_Seidel(v, niter, dx)

v = y3.copy()
res3Jacobi = Jacobi(v, niter, dx)

v = y3.copy()
res3GS = Gauss_Seidel(v, niter, dx)

v = y3.copy()
res3RBGS = Red_Black_Gauss_Seidel(v, niter, dx)

v = y6.copy()
res6Jacobi = Jacobi(v, niter, dx)

v = y6.copy()
res6GS = Gauss_Seidel(v, niter, dx)

v = y6.copy()
res6RBGS = Red_Black_Gauss_Seidel(v, niter, dx)

id = np.zeros(niter)
for i in range(0, niter):
   id[i] = i
   
#print("v=\n",v)   

plt.figure("OneFLOW-CFD Solver Jacobi+Gauss Seidel Residual", figsize=(6, 4), dpi=100)  
plt.plot(id,res1Jacobi,'r-',linewidth=2,label='k=1 Jacobi')
plt.plot(id,res1GS,'b-',linewidth=2,label='k=1 GS')
plt.plot(id,res1RBGS,'k-',linewidth=2,label='k=1 RBGS')
plt.plot(id,res3Jacobi,'r-.',linewidth=2,label='k=3 Jacobi')
plt.plot(id,res3GS,'b-.',linewidth=2,label='k=3 GS')
plt.plot(id,res3RBGS,'k-.',linewidth=2,label='k=3 RBGS')
plt.plot(id,res6Jacobi,'r--',linewidth=2,label='k=6 Jacobi')
plt.plot(id,res6GS,'b--',linewidth=2,label='k=6 GS')
plt.plot(id,res6RBGS,'k--',linewidth=2,label='k=6 RBGS')
plt.legend()
plt.show()
