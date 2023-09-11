import matplotlib.pyplot as plt
import numpy as np

n=8;
print("n=",n);
nx = n;
ny = n;

MaxStep=60;
Visc=0.1;
dt=0.02; 
# parametersforSOR iteration
MaxIt=100;
Beta=1.5;
MaxErr=0.001;
sf=np.zeros((nx+1,ny+1));
vt=np.zeros((nx+1,ny+1));
w=np.zeros((nx+1,ny+1));
h=1.0/(nx);
t=0.0;

uwall = 1.0
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
    x[i+1] = xmin + ( i + 1 ) * dx

for j in range(ny):
    y[j+1] = ymin + ( j + 1 ) * dy
    
  
# Creating 2-D grid of features
[X, Y] = np.meshgrid(x, y, indexing='ij')

fig = plt.figure("OneFLOW-CFD ", figsize=(10, 6), dpi=100)
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

h2 = h * h

print("h=",h,"h2=",h2)

for istep in range(MaxStep):
    print("istep=",istep+1)
    #solve for the streamfunction
    for iter in range(MaxIt):
        w = sf.copy(); # by SOR iteration
        for j in range(1,ny):
            for i in range(1,nx):
                sf[i,j] = 0.25 * Beta * ( sf[i+1,j] + sf[i-1,j] + sf[i,j+1] + sf[i,j-1] + h2 * vt[i,j] ) \
                        + ( 1.0 - Beta ) * sf[i,j];
        Err=0.0;
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                Err = Err + np.abs( w[i,j] - sf[i,j] );
        #stop if iteration has converged
        if Err<= MaxErr :
            break
            
    vt[1:nx,0 ] = -2.0 * sf[1:nx,   1]/h2 + 0 * 2.0/h;      # vorticity on bottom wall
    vt[1:nx,ny] = -2.0 * sf[1:nx,ny-1]/h2 - uwall * 2.0/h;  # vorticity on top wall
    vt[0 ,1:ny] = -2.0 * sf[1   ,1:ny]/h2;                  # vorticity on left wall
    vt[nx,1:ny] = -2.0 * sf[nx-1,1:ny]/h2;                  # vorticity on right wall

    for j in range(1,ny):
        for i in range(1,nx):
            w[i,j] = 0.25 * ( ( sf[i+1,j] - sf[i-1,j] ) * ( vt[i,j+1] - vt[i,j-1] ) - \
                              ( sf[i,j+1] - sf[i,j-1] ) * ( vt[i+1,j] - vt[i-1,j] ) ) / h2 \
                   + Visc * ( vt[i+1,j] + vt[i-1,j] + vt[i,j+1] + vt[i,j-1] - 4.0 * vt[i,j] ) / h2;
                   
    vt[1:nx,1:ny] = vt[1:nx,1:ny] + dt * w[1:nx,1:ny]; #update the vorticity

    t = t + dt
    #print("t=",t)
    
    ax1.clear()
    ax2.clear()
    ax1.contour(X, Y, vt)
    ax2.contour(X, Y, sf)
    plt.pause( 0.01 )


plt.show()