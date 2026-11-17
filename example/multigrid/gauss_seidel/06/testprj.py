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
    
mylabel_list=np.array(["Eigenvectors of RG Init","Eigenvectors of A Init"])
mycolor_list=np.array(['blue','red'])
mymarker_list=np.array(['o','s'])

nCase = len( mylabel_list )
print("nCase=",nCase)

iters = np.zeros( (nCase,n+1) )
wavenumbers = np.zeros( n+1 )

for k in range(1, n):
   wavenumbers[k] = k
   km = wavenumbers[k]
   aa = np.cos( km * np.pi / n )
   bb = 1.0
   for j in range(0, n+1):
       y[j] = bb * np.sin( j * km * np.pi / n )
       bb *= aa
   v = y.copy()
   iter = Gauss_Seidel_Iter_Number(v, 0.01, dx)
   if iter >= 100:
       iter = 100
   iters[0][k] = iter  

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
   iters[1][k] = iter  

fig = plt.figure("OneFLOW-CFD Solver Gauss–Seidel iteration VS Wavenumber k", figsize=(6, 8), dpi=100)

for j in range(0, nCase):
    ax = fig.add_subplot(nCase,1,j+1)
    ax.set_xlabel("Wavenumber k")
    ax.set_ylabel("Iterations")
    #ax.plot(wavenumbers[1:-2],iters[j][1:-2],color=mycolor_list[j],linewidth=2,marker=mymarker_list[j],markerfacecolor='none',label=mylabel_list[j])
    ax.plot(wavenumbers[1:-2],iters[j][1:-2],color=mycolor_list[j],linewidth=2,markerfacecolor='none',label=mylabel_list[j])
    ax.legend()

fig.tight_layout()

plt.show()
   
