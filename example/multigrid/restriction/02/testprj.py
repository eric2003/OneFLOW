import matplotlib.pyplot as plt
import numpy as np

def plot_point( xp, yp ):
    x = [xp]
    y = [yp]
    plt.plot(x, y, marker="o", markersize=8, markeredgecolor="blue",
    markerfacecolor="blue")
    
def FieldNorm(vh): 
    return np.max(np.abs(vh))
    
def Jacobi(v, f, n, niter, omega):
    v2 = v.copy()
    a = 0.5 * omega
    b = 1.0 - omega
    c = a    
    d = a
    for iter in range(0, niter):
        for i in range(1, n):
            v2[i] = a * v[i-1] + b * v[i] + c * v[i+1] + d * f[i]
        v = v2.copy()
    return v    
    
def gauss_seidel(un, f, r, nx, ny, dx, dy):
    dx2 = dx * dx
    dy2 = dy * dy
    den = -2.0/dx2 - 2.0/dy2
    for j in range(1, ny):
        for i in range(1, nx):
            d2udx2 = (un[i-1,j] - 2*un[i,j] + un[i+1,j])/(dx2)
            d2udy2 = (un[i,j-1] - 2*un[i,j] + un[i,j+1])/(dy2)
            r[i,j] = f[i,j] - d2udx2 - d2udy2
            un[i,j] = un[i,j] + r[i,j] / den
    

#n=4
n=16
nx = n
ny = n
h=1/n
print("n=",n)
print("h=",h)
dx = h
dy = h

x = np.zeros(n+1)
y = np.zeros(n+1)
xmin = 0
xmax = 1
ymin = 0
ymax = 1
x[0] = xmin
for i in range(n):
    #print("i=",i,"n=",n)
    x[i+1] = xmin + ( i + 1 ) * dx

for j in range(n):
    y[j+1] = ymin + ( j + 1 ) * dy
    
ue = np.zeros( ( n + 1, n + 1 ) )
un = np.zeros( ( n + 1, n + 1 ) )
f = np.zeros( ( n + 1, n + 1 ) )
r = np.zeros( ( n + 1, n + 1 ) )

for j in range(1, n):
    for i in range(0, n+1):
        y16m = np.sin( i * 16 * np.pi / n )
        y40m = np.sin( i * 40 * np.pi / n )
        #y16m = np.sin( i * np.pi / n )
        #y40m = np.sin( i * np.pi / n )
        ue[i,j] = 0.5 * ( y16m + y40m )
        un[i,j] = ue[i,j]
        
nCycle = 1000

res = np.zeros( nCycle + 1 )
res[0] = FieldNorm( un )

for iCycle in range(0, nCycle):
    gauss_seidel(un, f, r, nx, ny, dx, dy)
    res[iCycle+1] = FieldNorm( un )
    
    
print("res=\n",res)        
    
fig = plt.figure("OneFLOW-CFD Solver+Restriction", figsize=(6, 8), dpi=100)
plt.axis('off')

ax = fig.add_subplot(2, 2, 1)

plt.title( "Line Value" )

for j in range(n+1):
    print("j=",j,"un=\n",un[:,j])
    ax.plot(x,un[:,j],color='black',linewidth=1)

ax = fig.add_subplot(2, 2, 2)

for j in range(n+1):
    ax.hlines(y = y[j], xmin = x[0], xmax = x[-1], colors="k" )
    for i in range(n+1):
        text = "({0},{1})".format(i,j)
        #plt.text(x[i], y[j], text, fontsize=10)
        #plot_point( x[i], y[j] )

for i in range(n+1):
    ax.vlines(x = x[i], ymin = y[0], ymax = y[-1], colors="k")
    

plt.contourf(x, y, un, levels=20, cmap="jet")
#plt.contour(x, y, ue, levels=20, cmap="jet")
plt.colorbar()
plt.title( "Exact solution" )

ax = fig.add_subplot(2, 2, 3)
ii = np.linspace(0, nCycle, nCycle + 1)
ax.set_yscale("log", base=10)
ax.plot(ii,res,color='black',linewidth=2,linestyle='solid',marker='none',markerfacecolor='none',label='res')

plt.tight_layout()
plt.show()
