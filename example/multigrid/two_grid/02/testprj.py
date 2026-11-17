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
    return v

n = 64
m = int(n/2)
dx = 1.0 / n
h = dx

xh = np.zeros( n + 1 )
yh = np.zeros( n + 1 )
ih = np.zeros( n + 1 )

x = np.zeros( n + 1 )
y = np.zeros( n + 1 )

xh_odd = np.zeros( m + 1 )
yh_odd = np.zeros( m + 1 )
ih_odd = np.zeros( m + 1 )

xh_even = np.zeros( m + 1 )
yh_even = np.zeros( m + 1 )
ih_even = np.zeros( m + 1 )

x2h = np.zeros( m + 1 )
y2h = np.zeros( m + 1 )
i2h = np.zeros( m + 1 )

for i in range(0, n+1):
    x[i] = i * h
    y16m = np.sin( i * 16 * np.pi / n )
    y40m = np.sin( i * 40 * np.pi / n )
    y[i] = 0.5 * ( y16m + y40m )
    
omega = 2.0/3.0    
niter = 1
v = y.copy()
res = Weighted_Jacobi(v, niter, omega)
 
km = 2
io = 0
ie = 0
for j in range(0, n + 1):
    xh[j] = j * h
    ih[j] = j
    yh[j] = np.sin( j * km * np.pi / n )
    if (j % 2) == 0:
        xh_odd[io] = xh[j]
        ih_odd[io] = ih[j]
        yh_odd[io] = yh[j]
        io += 1
    else:
        xh_even[ie] = xh[j]
        ih_even[ie] = ih[j]
        yh_even[ie] = yh[j]
        ie += 1

#print("ih_odd=\n",ih_odd)
#print("ih_even=\n",ih_even)

for j in range(0, m + 1):
    i2h[j] = j
    x2h[j] = j * 2 * h
    y2h[j] = np.sin( j * km * np.pi / m )

fig = plt.figure('OneFLOW-CFD Solver Interpolation of a vector on coarse grid to fine grid ', figsize=(10, 4), dpi=100)

ax = fig.add_subplot(1,2,1)
#ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel("X")
ax.set_ylabel("Y",rotation='horizontal')
ax.xaxis.set_label_coords(1.01, 0.5)
ax.yaxis.set_label_coords(0.0, 1.0)
ax.set_ylim(-1.0, 1.0)

mycolor  = 'blue'
#mymarker = 'o'
mymarker = 'none'
mylabel  = 'coarse grid'
ax.plot(x,y,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(1,2,2)
#ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel("X")
ax.set_ylabel("Y",rotation='horizontal')
ax.xaxis.set_label_coords(1.01, 0.5)
ax.yaxis.set_label_coords(0.0, 1.0)
ax.set_ylim(-1.0, 1.0)

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
#ax.plot(x,res,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)
ax.plot(x,res,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

fig.tight_layout()
plt.show()
