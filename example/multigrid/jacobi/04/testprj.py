import numpy as np
import matplotlib.pyplot as plt

def Weighted_Jacobi(v, niter, omega):
    v2 = v.copy()
    coef = 1.0/3.0
    res = np.zeros(niter)
    a = 0.5 * omega
    c = a
    b = 1.0 - omega
    for iter in range(0, niter):
        for i in range(1, n):
            v2[i] = a * v[i-1] + b * v[i] + c * v[i+1]
        v = v2.copy()
        resm = np.abs(np.max(v2))
        res[iter] = resm
    return res

n = 64
h = 1.0/n
x = np.zeros( n+1 )
y1 = np.zeros( n+1 )
y3 = np.zeros( n+1 )
y6 = np.zeros( n+1 )
v = np.zeros( n+1 )

for i in range(0, n+1):
    x[i] = i * h
    y1[i] = np.sin( i * 1 * np.pi / n )
    y3[i] = np.sin( i * 3 * np.pi / n )
    y6[i] = np.sin( i * 6 * np.pi / n )    
    
niter = 100
print("v=\n",v)
omega = 2.0/3.0

v = y1.copy()
res1 = Weighted_Jacobi(v, niter, omega)

v = y3.copy()
res3 = Weighted_Jacobi(v, niter, omega)

v = y6.copy()
res6 = Weighted_Jacobi(v, niter, omega)

id = np.zeros(niter)
for i in range(0, niter):
   id[i] = i
   
fig = plt.figure("OneFLOW-CFD Solver Weighted Jacobi iteration Residual", figsize=(6, 4), dpi=100)

ax = fig.add_subplot()
ax.set_xlabel("Iterations")
ax.set_ylabel("error")

ax.plot(id,res1,'k-',linewidth=2,label='k=1')
ax.plot(id,res3,'b-.',linewidth=2,label='k=3')
ax.plot(id,res6,'g--',linewidth=2,label='k=6')
ax.legend()
plt.show()



