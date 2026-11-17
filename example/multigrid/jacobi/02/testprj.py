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
y3 = np.zeros( n+1 )
y6 = np.zeros( n+1 )
v = np.zeros( n+1 )

for i in range(0, n+1):
    x[i] = i * h
    y1[i] = np.sin( i * 1 * np.pi / n )
    y3[i] = np.sin( i * 3 * np.pi / n )
    y6[i] = np.sin( i * 6 * np.pi / n )
    
#print("x=\n",x)


#print("v=\n",v)
#resv = np.abs(np.max(v))
#print("resv=",resv)

niter = 100

v = y1.copy()
res1 = Jacob(v, niter)

v = y3.copy()
res3 = Jacob(v, niter)

v = y6.copy()
res6 = Jacob(v, niter)

id = np.zeros(niter)
for i in range(0, niter):
   id[i] = i
   
#print("res=\n",res)   

fig = plt.figure("OneFLOW-CFD Solver Weighted Jacobi iteration Residual", figsize=(6, 4), dpi=100)

ax = fig.add_subplot()
ax.set_yscale("log", base=10)
ax.set_xlabel("Iterations")
ax.set_ylabel("Log error")
ax.set_ylim([1e-3, 2.0])

ax.plot(id,res1,'k-',linewidth=2,label='k=1')
ax.plot(id,res3,'b-.',linewidth=2,label='k=3')
ax.plot(id,res6,'g--',linewidth=2,label='k=6')
ax.legend()
plt.show()


