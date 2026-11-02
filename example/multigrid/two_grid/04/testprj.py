import numpy as np
import matplotlib.pyplot as plt

def Weighted_Jacobi(v, f, n, niter, omega):
    v2 = v.copy()
    res = np.zeros(niter)
    a = 0.5 * omega
    b = 1.0 - omega
    c = a    
    d = a
    for iter in range(0, niter):
        for i in range(1, n):
            v2[i] = a * v[i-1] + b * v[i] + c * v[i+1] + d * f[i]
        v = v2.copy()
    return v

def CalcResidual(v, f, n):
    v2 = v.copy()
    a = -1.0
    b = 2.0
    c = -1.0
    for i in range(1, n):
        av = a * v[i-1] +  b * v[i] + c * v[i+1]
        v2[i] = f[i] - av
    return v2
    
def Restrict(rh, nc):
    r2h = np.zeros( nc + 1 )
    for i in range(1, nc):
        r2h[i] = rh[2*i]
    return r2h

def Interpolate(e2h, nc):
    eh = np.zeros( 2*nc + 1 )
    for i in range(1, nc):
        eh[2*i] = e2h[i]
        eh[2*i+1] = 0.5 * ( e2h[i] + e2h[i+1] )
    return eh
    
def Update(vh, eh, n):
    vh_new = np.zeros( n + 1 )
    for i in range(1, n):
        vh_new[i] = vh[i] + eh[i]
    return vh_new

def basic_visual_setup(ax):    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlabel("X")
    ax.set_ylabel("Y",rotation='horizontal')
    ax.xaxis.set_label_coords(1.02, 0.02)
    ax.yaxis.set_label_coords(0.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    return


n = 64
m = int(n/2)
dx = 1.0 / n
h = dx

xh = np.zeros( n + 1 )
yh = np.zeros( n + 1 )
ih = np.zeros( n + 1 )
fh = np.zeros( n + 1 )
rh = np.zeros( n + 1 )
eh = np.zeros( n + 1 )

x = np.zeros( n + 1 )
y = np.zeros( n + 1 )

x2h = np.zeros( m + 1 )
y2h = np.zeros( m + 1 )
e2h = np.zeros( m + 1 )

for i in range(0, n+1):
    x[i] = i * h
    y16m = np.sin( i * 16 * np.pi / n )
    y40m = np.sin( i * 40 * np.pi / n )
    y[i] = 0.5 * ( y16m + y40m )
    
omega = 2.0/3.0
vh = y.copy()
vh_1 = Weighted_Jacobi(vh, fh, n, 1, omega)
vh_3 = Weighted_Jacobi(vh, fh, n, 3, omega)
rh  = CalcResidual(vh_3, fh, n)
r2h = Restrict(rh, int(n/2) )
x2h = Restrict(x, int(n/2) )

e2h = Weighted_Jacobi(e2h, r2h, int(n/2), 3, omega)
eh  = Interpolate(e2h, int(n/2) )

vh_new = Update(vh_3, eh, n)
vh_4 = Weighted_Jacobi(vh_new, fh, n, 3, omega)

fig = plt.figure('OneFLOW-CFD Solver Two-Grid Correction Scheme', figsize=(10, 8), dpi=100)

ax = fig.add_subplot(3,2,1)
basic_visual_setup( ax )
mycolor  = 'blue'
mymarker = 'none'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(3,2,2)
basic_visual_setup( ax )

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
ax.plot(x,vh_1,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(3,2,3)
basic_visual_setup( ax )

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
ax.plot(x,vh_3,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(3,2,4)
basic_visual_setup( ax )

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
ax.plot(x,vh_new,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(3,2,5)
basic_visual_setup( ax )

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
ax.plot(x,eh,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(3,2,6)
basic_visual_setup( ax )

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
ax.plot(x,vh_4,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

fig.tight_layout()
plt.show()
