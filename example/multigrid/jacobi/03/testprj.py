import numpy as np
import matplotlib.pyplot as plt

def Jacob(v, niter):
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

n = 64
h = 1.0/n
x = np.zeros( n+1 )
y1 = np.zeros( n+1 )
y6 = np.zeros( n+1 )
y32 = np.zeros( n+1 )
v = np.zeros( n+1 )

for i in range(0, n+1):
    x[i] = i * h
    y1[i] = np.sin( i * 1 * np.pi / n )
    y6[i] = np.sin( i * 6 * np.pi / n )
    y32[i] = np.sin( i * 32 * np.pi / n )
    v[i] = (y1[i]+y6[i]+y32[i])/3.0
    
niter = 100
print("v=\n",v)
res = Jacob(v, niter)

id = np.zeros(niter)
for i in range(0, niter):
   id[i] = i
   
print("res=\n",res)   

fig = plt.figure("OneFLOW-CFD Solver Weighted Jacobi iteration Residual", figsize=(6, 4), dpi=100)

ax = fig.add_subplot()
ax.set_xlabel("Iterations")
ax.set_ylabel("error")

ax.plot(id,res,'b-',linewidth=2,label='k=(v1+v6+v32)/3')
ax.legend()
plt.show()


