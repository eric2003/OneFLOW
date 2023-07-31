import numpy as np
import matplotlib.pyplot as plt

def Weighted_Jacobi(v, niter, omega):
    v2 = v.copy()
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
    
def Weighted_Jacobi_Iter_Number(v, factor, omega):
    v2 = v.copy()
    a = 0.5 * omega
    c = a
    b = 1.0 - omega
    res0 = np.abs(np.max(v2))
    res_ratio = 1.0
    iter_number = 0
    while res_ratio > factor:
        for i in range(1, n):
            v2[i] = a * v[i-1] + b * v[i] + c * v[i+1]
        v = v2.copy()
        resm = np.abs(np.max(v2))
        res_ratio = resm / res0
        print("iter=,ratio=",iter_number,res_ratio)
        iter_number += 1
    return iter_number    

n = 64
h = 1.0/n
x = np.zeros( n+1 )
y = np.zeros( n+1 )
v = np.zeros( n+1 )


omega = 1.0
mylabel='w=1'
#omega = 2.0/3.0
#mylabel='w=2/3'
iters = np.zeros( n+1 )
wavenumbers = np.zeros( n+1 )

for k in range(1, n):
    wavenumbers[k] = k
    km = wavenumbers[k]
    for i in range(0, n+1):
        x[i] = i * h
        y[i] = np.sin( i * km * np.pi / n )
    v = y.copy()
    iter = Weighted_Jacobi_Iter_Number(v, 0.01, omega)
    if iter >= 100:
        iter = 100
    iters[k] = iter
    
fig = plt.figure("OneFLOW-CFD Solver Weighted Jacobi iteration VS Wavenumber k", figsize=(6, 4), dpi=100)

ax = fig.add_subplot()
ax.set_xlabel("Wavenumber k")
ax.set_ylabel("Iterations")

ax.plot(wavenumbers[1:-2],iters[1:-2],'b-',linewidth=2,label=mylabel)
ax.legend()
plt.show()
   

