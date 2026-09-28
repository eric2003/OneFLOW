import numpy as np
import matplotlib.pyplot as plt
import time

def compute_l2norm( nx, ny, r ):
    rms = 0.0
    for j in range(0, ny+1):
        for i in range(0, nx+1):
            rms = rms + r[i,j] * r[i,j]
    rms = np.sqrt( rms / ( ( nx + 1 ) * ( ny + 1 ) ) )
    return rms


def my_fst( nx, ny, f, fmn, aim, bjn ):
    for n in range(0, ny):
        print("my_fst n,ny=",n,ny)    
        for m in range(0, nx):
            fmn[m,n] = 0.0
            for j in range(0, ny):
                c1 = bjn[j,n]
                for i in range(0, nx):
                    c2 = aim[i,m]
                    fmn[m,n] += f[i,j] * c1 * c2
    return

def inverse_fst( nx, ny, umn, aim, bjn ):
    uij = np.zeros( ( nx, ny ) )
    coef = 4.0 / ( nx * ny )
    for j in range(0, ny):
        print("inverse_fst j=,ny=",i,ny)
        for i in range(0, nx):    
            uij[i,j] = 0.0
            for n in range(0, ny):
                c1 = bjn[j,n]
                for m in range(0, nx):            
                    c2 = aim[i,m]
                    uij[i,j] += umn[m,n] * c1 * c2
            uij[i,j] *= coef
    return uij

def compute_matrix( nx, ny, aim, bjn ):
    for i in range(0, nx):    
        for m in range(0, nx):
            aim[i,m] = np.sin(np.pi*m*i/nx)

    for j in range(0, ny):  
        for n in range(0, ny):
            bjn[j,n] = np.sin(np.pi*n*j/ny)
    return

def poisson_fst( nx, ny, dx, dy, f ):
    uij = np.zeros( ( nx, ny ) )
    aim = np.zeros( ( nx, ny ) )
    bjn = np.zeros( ( nx, ny ) )
    compute_matrix( nx, ny, aim, bjn )
    aa = -2.0/(dx*dx) - 2.0/(dy*dy)
    bb = 2.0/(dx*dx)
    cc = 2.0/(dy*dy)
    
    fmn = np.zeros( ( nx, ny ) )
    umn = np.zeros( ( nx, ny ) )
    my_fst( nx, ny, f, fmn, aim, bjn )
    
    small = 1.0e-10

    for j in range(0, ny):
        for i in range(0, nx):
            umn[i,j] = fmn[i,j] / ( aa + bb*np.cos(np.pi * i / nx) + cc*np.cos(np.pi * j / ny) + small )

    uij = inverse_fst( nx, ny, umn, aim, bjn )

    return uij

nx = 512
ny = 512

#nx = 128
#ny = 128

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

dx = (x_r-x_l)/nx
dy = (y_t-y_b)/ny

x = np.zeros( nx + 1 )
y = np.zeros( ny + 1 )
ue = np.zeros( ( nx + 1, ny + 1 ) )
f  = np.zeros( ( nx + 1, ny + 1 ) )
un = np.zeros( ( nx + 1, ny + 1 ) )

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)

for j in range(0, ny+1):
    y[j] = y_b + dy*(j)

# given exact solution
km = 16
c1 = (1.0/km)**2
c2 = -8.0*np.pi*np.pi 

for j in range(0, ny+1):
    for i in range(0, nx+1):
        ue[i,j] = np.sin(2.0*np.pi*x[i])*np.sin(2.0*np.pi*y[j]) + \
                  c1*np.sin(km*2.0*np.pi*x[i])*np.sin(km*2.0*np.pi*y[j])
               
        f[i,j] = c2*np.sin(2.0*np.pi*x[i])*np.sin(2.0*np.pi*y[j]) + \
                 c2*np.sin(km*2.0*np.pi*x[i])*np.sin(km*2.0*np.pi*y[j])

        un[i,j] = 0.0
        
#print("nx=",nx)
#print("ny=",ny)
#print("len(ue[0,:]=",len(ue[0,:]))
#print(ue)        
        
#plt.figure("OneFLOW-CFD Solver+2D Poisson Equation+FST", figsize=(5, 4), dpi=100)
#plt.imshow(ue, extent=[0, 1, 0, 1], cmap="jet", origin="lower")
#plt.colorbar()
#plt.tight_layout() 
#plt.show()       
        
time_1 = time.time()
tmp = poisson_fst(nx,ny,dx,dy,f)
for j in range(0, ny):
    for i in range(0, nx):
        un[i,j] = tmp[i,j]
# Periodic boundary condition
un[nx,:] = un[0,:]
un[:,ny] = un[:,0]

time_2 = time.time()
t = time_2 - time_1

uerror = np.zeros( ( nx + 1, ny + 1 ) )
rms_error = 0.0

uerror = un - ue

rms_error = compute_l2norm(nx, ny, uerror)
max_error = np.max( np.abs(uerror) )

print("Error details:");
print("L-2 Norm = ", rms_error);
print("Maximum Norm = ", max_error);
print("CPU Time = ", t);

output = open("output.txt", "w");
output.write("Error details: \n");
output.write("L-2 Norm = {0}\n" .format(rms_error));
output.write("Maximum Norm = {0}\n".format(max_error));
output.write("CPU Time = {0}\n".format(t));
output.close()

titles = ("Exact solution", "Numerical solution")

plt.figure("OneFLOW-CFD Solver+2D Poisson Equation+FST", figsize=(10, 4), dpi=100)
for i, result in enumerate([ue, un]):
    plt.subplot(1, 2, i+1)
    plt.imshow(result, extent=[0, 1, 0, 1], cmap="jet", origin="lower")
    plt.colorbar()
    plt.title(titles[i])
plt.tight_layout()
plt.show()
