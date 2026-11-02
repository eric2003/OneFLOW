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

    
def ps_fft( nx, ny, dx, dy, f ):
    eps = 1.0e-6

    kx = np.zeros( nx )
    ky = np.zeros( ny )

    data  = np.zeros( ( nx, ny ), dtype=complex )
    data1 = np.zeros( ( nx, ny ), dtype=complex )
    e = np.zeros( ( nx, ny ) )
    u = np.zeros( ( nx, ny ) )

    aa = -2.0/(dx*dx) - 2.0/(dy*dy)
    bb = 2.0/(dx*dx)
    cc = 2.0/(dy*dy)

    #wave number indexing
    hx = 2 * np.pi / nx
    nn = int(nx / 2)
    for i in range(0,nn):
        kx[i] = hx * i
        kx[i + nn] = hx * ( i - nn )        

    kx[0] = eps
    ky = kx.copy()

    for j in range(0, ny):
        for i in range(0, nx):
            data[i,j] = complex(f[i,j],0.0)

    e = np.fft.fft2(data)
    e[0,0] = 0.0
    
    for j in range(0, ny):
        for i in range(0, nx):
            data1[i,j] = e[i,j] / ( aa + bb*np.cos(kx[i]) + cc*np.cos(ky[j]) )

    u = np.fft.ifft2(data1).real

    return u
    
def fft( nx, ny ):
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
            
    tmp = ps_fft(nx,ny,dx,dy,f)
    for j in range(0, ny):
        for i in range(0, nx):
            un[i,j] = tmp[i,j]
    # Periodic boundary condition
    un[nx,:] = un[0,:]
    un[:,ny] = un[:,0]
    
    uerror = np.zeros( ( nx + 1, ny + 1 ) )
    rms_error = 0.0
    
    uerror = un - ue
    
    rms_error = compute_l2norm(nx, ny, uerror)
    max_error = np.max( np.abs(uerror) )
    
    print("Error details:");
    print("L-2 Norm = ", rms_error);
    print("Maximum Norm = ", max_error);
    return rms_error    

res = open("center_ftt_res.txt", "w");

grid = [32,64,128,256,512]
for i in range( 0, len( grid ) ):
    nx = grid[i]
    ny = grid[i]
    error = fft( nx, ny )
    print("i=",i, " nx = ", nx, " ny = ", ny, "error=", error)
    res.write("{0} {1}\n" .format(grid[i], error));
res.close()

