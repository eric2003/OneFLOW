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
        print("iter=",iter,"resm=",resm)
        res[iter] = resm
    return v

n = 64
h = 1.0/n
x = np.zeros( n+1 )

omega = 2.0/3.0
mylabel='w=2/3'
iters = np.zeros( n+1 )
wavenumbers = np.zeros( n+1 )

y3 = np.zeros( n+1 )
y16 = np.zeros( n+1 )
y2_16 = np.zeros( n+1 )

yy = np.zeros( (6,n+1) )

for i in range(0, n+1):
    x[i] = i * h
    y2m = np.sin( i * 2 * np.pi / n )
    y3m = np.sin( i * 3 * np.pi / n )
    y16m = np.sin( i * 16 * np.pi / n )
    y3[i] = y3m
    y16[i] = y16m
    y2_16[i] = 0.5 * ( y2m + y16m )
    
yy[0] = y3.copy()
yy[1] = y3.copy()

yy[2] = y16.copy()
yy[3] = y16.copy()

yy[4] = y2_16.copy()
yy[5] = y2_16.copy()

iters= np.array([1,10,1,10,1,10])

for k in range(0, 6):
    yy[k] = Weighted_Jacobi(yy[k],iters[k], omega)

fig = plt.figure("OneFLOW-CFD Solver Weighted Jacobi iteration with different wavenumber", figsize=(8, 6), dpi=100)

for k in range(0, 6):
    ax = fig.add_subplot(3,2,k+1)
    ax.set_xlabel("x")
    ax.set_ylabel("v")
    ax.set_ylim([-1.05, 1.05])
    ax.plot(x,yy[k],'k-',linewidth=2,label=mylabel)
    #ax.legend()

fig.tight_layout()
plt.show()

