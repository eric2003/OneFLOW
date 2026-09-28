import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
def compute_l2norm( nx, ny, r ):
    rms = 0.0
    for j in range(0, ny+1):
        for i in range(0, nx+1):
            rms = rms + r[i,j] * r[i,j]
    rms = np.sqrt( rms / ( ( nx + 1 ) * ( ny + 1 ) ) )
    return rms
    
#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -J(w,ψ) + ν ∇^2(w)
#-----------------------------------------------------------------------------#
def rhs(nx,ny,dx,dy,re,w,s,r):
    # Arakawa numerical scheme for Jacobian
    aa = 1.0/(re*dx*dx)
    bb = 1.0/(re*dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0
    
    for j in range(1, ny):
        for i in range(1, nx):
            j1 = gg*((w[i+1,j]-w[i-1,j])*(s[i,j+1]-s[i,j-1]) - \
                     (w[i,j+1]-w[i,j-1])*(s[i+1,j]-s[i-1,j]))
            
            j2 = gg*(w[i+1,j]*(s[i+1,j+1]-s[i+1,j-1]) - \
                     w[i-1,j]*(s[i-1,j+1]-s[i-1,j-1]) - \
                     w[i,j+1]*(s[i+1,j+1]-s[i-1,j+1]) + \
                     w[i,j-1]*(s[i+1,j-1]-s[i-1,j-1]))
            
            j3 = gg*(w[i+1,j+1]*(s[i,j+1]-s[i+1,j]) - \
                     w[i-1,j-1]*(s[i-1,j]-s[i,j-1]) - \
                     w[i-1,j+1]*(s[i,j+1]-s[i-1,j]) + \
                     w[i+1,j-1]*(s[i+1,j]-s[i,j-1]))
            
            jac = (j1+j2+j3)*hh
            
            #Central difference for Laplacian
            r[i,j] = -jac + aa*(w[i+1,j]-2.0*w[i,j]+w[i-1,j]) + \
                            bb*(w[i,j+1]-2.0*w[i,j]+w[i,j-1])
                            
                            
def bc(nx,ny,dx,dy,w,s):
    # first order approximation
    # boundary condition for vorticity (Hoffmann) left and right
    for j in range(0, ny+1):
        w[0 ,j] = -2.0*s[1   ,j]/(dx*dx)
        w[nx,j] = -2.0*s[nx-1,j]/(dx*dx)

    # boundary condition for vorticity (Hoffmann) bottom and top
    for i in range(0, nx+1):
        w[i,0 ] = -2.0*s[i,1   ]/(dy*dy)
        w[i,ny] = -2.0*s[i,ny-1]/(dy*dy) - 2.0/dy
                        
def bc2(nx,ny,dx,dy,w,s):
    # second order approximation
    # boundary condition for vorticity (Jensen) left and right
    for j in range(0, ny+1):
        w[0 ,j] = (-4.0*s[1   ,j]+0.5*s[2   ,j])/(dx*dx)
        w[nx,j] = (-4.0*s[nx-1,j]+0.5*s[nx-2,j])/(dx*dx)

    # boundary condition for vorticity (Jensen) bottom and top
    for i in range(0, nx+1):
        w[i,0 ] = (-4.0*s[i,1   ]+0.5*s[i,2   ])/(dy*dy)
        w[i,ny] = (-4.0*s[i,ny-1]+0.5*s[i,ny-2])/(dy*dy) - 3.0/dy
        
def cal_stream_function(nx,ny,dx,dy,sf,vt):
    # parameters for SOR iteration
    MaxIt=100;
    beta=1.5;
    MaxErr=1.0e-6;
    
    b = dx / dy
    b2 = b * b
    dx2 = dx * dx
    cc = 1.0 /( 2.0 * ( 1.0 + b2 ) )
    err0 = 1.0
    for iter in range(MaxIt):
        w = sf.copy(); # by SOR iteration
        for j in range(1,ny):
            for i in range(1,nx):
                sfnew = cc * ( dx2 * vt[i,j] + sf[i+1,j] + sf[i-1,j] + b2 * ( sf[i,j+1] + sf[i,j-1] ) )
                sf[i,j] = beta * sfnew + ( 1.0 - beta ) * sf[i,j];
        err=0.0;
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                err = err + np.abs( w[i,j] - sf[i,j] );
        if iter == 0:
           err0 = err
           #print("iter=",iter," err0=",err0)
           if err0 == 0:
              err0 = 1.0
        rerr = err / err0
        #print("iter=",iter," a=",a)
        #stop if iteration has converged
        if rerr <= 0.001 or err < MaxErr:
            break
                        
#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
def numerical(nx,ny,nt,ns,dx,dy,dt,re,wn):
    wt = np.zeros((nx+2,ny+2)) # temporary array during RK3 integration
    r  = np.zeros((nx+1,ny+1)) # right hand side
    sp = np.zeros((nx+1,ny+1)) # old streamfunction
    sn = np.zeros((nx+1,ny+1))
    ut = np.zeros((nx+1,ny+1)) 
    
    m = 1 # record index
    freq = int(nt/ns)
    print("freq=",freq)
    
    for k in range(0, nt):
        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,sn,r)
        
        for j in range(1, ny+1):
            for i in range(1, nx+1):
                wt[i,j] = wn[i,j] + dt*r[i,j]
                
        boundary(nx,ny,wt)

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for j in range(1, ny+1):
            for i in range(1, nx+1):
                wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]

        boundary(nx,ny,wt)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)
        
        for j in range(1, ny+1):
            for i in range(1, nx+1):
                wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]
                
        boundary(nx,ny,wn)
        
        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt)
        
        if ( k%freq == 0 ):
            print(k)
           
            for j in range(0, ny+1):
                for i in range(0, nx+1):   
                    ut[i,j] = wn[i+1,j+1]
                    
            filename = "vm"+str(m)+".txt"
            print("filename=",filename)
                    
            dump_field(x,y,ut,nx,ny,filename)
            
            m = m + 1


def boundary(nx,ny,f):
    # periodic BC
    for j in range(0, ny+2):
        f[0   , j] = f[nx, j]
        f[nx+1, j] = f[1 , j]

    for i in range(0, nx+2):
        f[i, 0   ] = f[i, ny]
        f[i, ny+1] = f[i, 1 ]
        
# initial condition for vortex merger problem
def vm_ic(nx,ny,x,y,w):
    sigma = np.pi
    xc1 = np.pi - np.pi/4.0
    yc1 = np.pi
    xc2 = np.pi + np.pi/4.0
    yc2 = np.pi
    
    for j in range(1, ny+1):
        for i in range(1, nx+1):    
            w[i,j] = np.exp(-sigma*((x[i-1]-xc1)**2 + (y[j-1]-yc1)**2)) + \
                     np.exp(-sigma*((x[i-1]-xc2)**2 + (y[j-1]-yc2)**2))
                     
def dump_field(x,y,field,nx,ny,filename):
    file = open(filename, "w")
    file.write("{0} {1}\n" .format(nx,ny));
    for j in range(0, ny+1):
        for i in range(0, nx+1):
            file.write("{0} {1} {2}\n" .format(x[i],y[j],field[i,j]));
    file.close()
                     

nx = 128
ny = 128

dt = 0.01
tf = 20.0
nt = int(tf/dt)
re = 1000.0
ns = 10

# parameters for SOR iteration
MaxIt=100;
Beta=1.5;
MaxErr=1.0e-6;

h=1.0/(nx);
t=0.0;

uwall = 1.0
x = np.zeros(nx+1)
y = np.zeros(ny+1)
rms = np.zeros(nt)

xmin = 0
xmax = 2.0 * np.pi
ymin = 0
ymax = 2.0 * np.pi

dx = ( xmax - xmin ) / nx
dy = ( ymax - ymin ) / ny
x[0] = xmin
y[0] = ymin
for i in range(nx):
    x[i+1] = xmin + ( i + 1 ) * dx

for j in range(ny):
    y[j+1] = ymin + ( j + 1 ) * dy

wn = np.zeros((nx+2,ny+2));
un = np.zeros((nx+1,ny+1));
un0 = np.zeros((nx+1,ny+1));
ue = np.zeros((nx+1,ny+1));
uerror = np.zeros((nx+1,ny+1));

time = 0.0
vm_ic(nx,ny,x,y,wn)

boundary(nx,ny,wn)
for j in range(0, ny+1):
    for i in range(0, nx+1):
        un0[i,j] = wn[i+1,j+1]
        
        
dump_field(x,y,un0,nx,ny,"vm0.txt")

numerical(nx,ny,nt,ns,dx,dy,dt,re,wn)

for j in range(0, ny+1):
    for i in range(0, nx+1):   
        un[i,j] = wn[i+1,j+1]

dump_field(x,y,un,nx,ny,"field_final.txt")
