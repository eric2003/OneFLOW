import matplotlib.pyplot as plt
import numpy as np

def plot_point( xp, yp, ax ):
    x = [xp]
    y = [yp]
    ax.plot(x, y, marker="o", markersize=8, markeredgecolor="blue",
    markerfacecolor="blue")

"""    
def restriction(nxf, nyf, nxc, nyc, r, ec):
    for j in range(1, nyc):
        for i in range(1, nxc):
            # grid index for fine grid for the same coarse point
            center = 4.0*r[2*i-1, 2*j-1]
            # E, W, N, S with respect to coarse grid point in fine grid
            grid = 2.0*(r[2*i-1, 2*j-1+1] + r[2*i-1, 2*j-1-1] +
                        r[2*i-1+1, 2*j-1] + r[2*i-1-1, 2*j-1])
            # NE, NW, SE, SW with respect to coarse grid point in fine grid
            corner = 1.0*(r[2*i-1+1, 2*j-1+1] + r[2*i-1+1, 2*j-1-1] +
                          r[2*i-1-1, 2*j-1+1] + r[2*i-1-1, 2*j-1-1])
            # restriction using trapezoidal rule
            ec[i,j] = (center + grid + corner)/16.0    
"""

def restriction(nxf, nyf, nxc, nyc):
    for j in range(1, nyc):
        for i in range(1, nxc):
            print("restriction i,j=",i,j)
            print("restriction 2*i-1,2*j-1=",2*i-1,2*j-1)


#n=128
n=4
print("n=",n)
nx = n
ny = n

x = np.zeros(nx+1)
y = np.zeros(ny+1)

xmin = 0
xmax = 1
ymin = 0
ymax = 1
dx = ( xmax - xmin ) / nx
dy = ( ymax - ymin ) / ny
x[0] = xmin
y[0] = ymin
for i in range(nx):
    print("i=",i,"nx=",nx)
    x[i+1] = xmin + ( i + 1 ) * dx

for j in range(ny):
    y[j+1] = ymin + ( j + 1 ) * dy
    
nxc = int(nx / 2)
nyc = int(ny / 2)
print("nxc=",nxc,"nyc=",nyc)
xc = np.zeros(nxc+1)
for ic in range(nxc+1):
    print("ic=",ic,"nxc=",nxc,"nx=",nx,"2*nxc",2*nxc)
    xc[ic] = x[2*ic]
    #xc[ic] = 0.5 * ( x[2*ic] + x[2*ic+1] )
    
print("xc=",xc)
yc = np.zeros(nyc+1)
for jc in range(nyc+1):
    yc[jc] = y[2*jc]
    
    
restriction(nx, ny, nxc, nyc)

fig = plt.figure("OneFLOW-CFD Solver+Restriction", figsize=(6, 8), dpi=100)

ax = fig.add_subplot(2, 1, 1)
ax.axis('off')

for j in range(ny+1):
    ax.hlines(y = y[j], xmin = x[0], xmax = x[-1], colors="k" )
    for i in range(nx+1):
        text = "({0},{1})".format(i,j)
        ax.text(x[i], y[j], text, fontsize=10)
        plot_point( x[i], y[j], ax )

for i in range(nx+1):
    ax.vlines(x = x[i], ymin = y[0], ymax = y[-1], colors="k")
    
ax = fig.add_subplot(2, 1, 2)
ax.axis('off')

for j in range(nyc+1):
    ax.hlines(y = yc[j], xmin = xc[0], xmax = xc[-1], colors="k" )
    for i in range(nxc+1):
        text = "({0},{1})".format(i,j)
        ax.text(xc[i], yc[j], text, fontsize=10)
        plot_point( xc[i], yc[j], ax )

for i in range(nxc+1):
    ax.vlines(x = xc[i], ymin = yc[0], ymax = yc[-1], colors="k")

plt.tight_layout()
plt.show()
