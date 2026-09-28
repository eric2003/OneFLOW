import numpy as np
import matplotlib.pyplot as plt

def Gauss_Seidel_Iter_Number(v, factor, dx):
    dx2 = dx * dx
    den = - 2.0 / dx2
    res0 = np.max(np.abs(v))
    #print( "res0=",res0)
    #print( "v=",v)
    res_ratio = 1.0
    iter_number = 0
    while res_ratio > factor:
        for i in range(1, n):
            d2udx2 = ( v[i-1] - 2*v[i] + v[i+1] )/dx2
            rm = -d2udx2
            v[i] = v[i] + rm / den
        resm = np.max(np.abs(v))
        res_ratio = resm / res0
        print("iter=,ratio=",iter_number,res_ratio)
        iter_number += 1
    return iter_number    
    
n = 64
h = 1.0/n
dx = h

print("dx=",dx)

x = np.zeros( n+1 )
y = np.zeros( n+1 )

for i in range(0, n+1):
    x[i] = i * dx
    
mylabel='w=2/3'
iters = np.zeros( n+1 )
wavenumbers = np.zeros( n+1 )

for k in range(1, n):
    wavenumbers[k] = k
    km = wavenumbers[k]
    aa = np.cos( km * np.pi / n )
    bb = 1.0
    for j in range(0, n+1):
        y[j] = bb * np.sin( j * km * np.pi / n )
        #bb *= aa
    v = y.copy()
    iter = Gauss_Seidel_Iter_Number(v, 0.01, dx)
    if iter >= 100:
        iter = 100
    iters[k] = iter  

fig = plt.figure("OneFLOW-CFD Solver Gauss–Seidel iteration VS Wavenumber k", figsize=(6, 4), dpi=100)

ax = fig.add_subplot()
ax.set_xlabel("Wavenumber k")
ax.set_ylabel("Iterations")

ax.plot(wavenumbers[1:-2],iters[1:-2],'b-',linewidth=2,label=mylabel)
ax.legend()
plt.show()
   
