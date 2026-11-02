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
                            
def rhsDEBUG(nx,ny,dx,dy,re,w,s,r):
    # Arakawa numerical scheme for Jacobian
    aa = 1.0/(re*dx*dx)
    bb = 1.0/(re*dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0
    
    for j in range(1, ny):
        for i in range(1, nx):
            j01 = (w[i+1,j]-w[i-1,j])*(s[i,j+1]-s[i,j-1])
            j02 = (w[i,j+1]-w[i,j-1])*(s[i+1,j]-s[i-1,j])
            j03 = (w[i,j+1]-w[i,j-1])
            j04 = (s[i+1,j]-s[i-1,j])
            j05 = s[i+1,j]
            j06 = s[i-1,j]
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
            c1 = aa*(w[i+1,j]-2.0*w[i,j]+w[i-1,j])
            c2 = bb*(w[i,j+1]-2.0*w[i,j]+w[i,j-1])
            if i==3 and j==3 :
                print("jac=",jac," c1=",c1, " c2=",c2)
                print("j1=",j1," j2=",j2, " j3=",j3)
                print("j01=",j01," j02=",j02)
                print("j03=",j03," j04=",j04)
                print("j05=",j05," j06=",j06)
            
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
        
def cal_stream_function(nx,ny,dx,dy,sf,vt,MaxIt,MaxErr,beta):
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
            
def cal_stream_function_DEBUG(nx,ny,dx,dy,sf,vt,MaxIt,MaxErr,beta):
    b = dx / dy
    b2 = b * b
    dx2 = dx * dx
    cc = 1.0 /( 2.0 * ( 1.0 + b2 ) )
    err0 = 1.0
    
    #for j in range(0,ny+1):
    #    for i in range(0,nx+1):
    #        print("i,j=",i, " ", j, " sf[i,j]=",sf[i,j])
    for iter in range(MaxIt):
        w = sf.copy(); # by SOR iteration
        for j in range(1,ny):
            for i in range(1,nx):
                sfnew = cc * ( dx2 * vt[i,j] + sf[i+1,j] + sf[i-1,j] + b2 * ( sf[i,j+1] + sf[i,j-1] ) )
                sf[i,j] = beta * sfnew + ( 1.0 - beta ) * sf[i,j];
                #print("i,j=",i, " ", j, " sf[i,j]=",sf[i,j])
        print("iter=",iter,"sf[2,3]=",sf[2,3])
        #exit();
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
        print("iter=",iter," err=",err, "MaxErr=", MaxErr)
        #stop if iteration has converged
        if rerr <= 0.001 or err < MaxErr:
            break            
                        
#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
def numericalold(nx,ny,nt,dx,dy,dt,re,wn,sn,rms,MaxIt,MaxErr,beta):
    wt = np.zeros((nx+1,ny+1)) # temporary array during RK3 integration
    r  = np.zeros((nx+1,ny+1)) # right hand side
    sp = np.zeros((nx+1,ny+1)) # old streamfunction
    
    for k in range(0, nt):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                sp[i,j] = sn[i,j]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,sn,r)
        
        for j in range(1, ny):
            for i in range(1, nx):
                wt[i,j] = wn[i,j] + dt*r[i,j]

        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for j in range(1, ny):
            for i in range(1, nx):
                wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]

        bc2(nx,ny,dx,dy,wt,sn)
        """
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " wt[i,j]=",wt[i,j])    
                
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " wn[i,j]=",wn[i,j])                    
                
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " r[i,j]=",r[i,j])                       
        exit()                
        """

        # compute streamfunction from vorticity
        #print(" sn[2,3]=",sn[2,3])    
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)
        #cal_stream_function_DEBUG(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        #rhsDEBUG(nx,ny,dx,dy,re,wt,sn,r)
        """
        print(" sn[2,3]=",sn[2,3]) 
        
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " r[i,j]=",r[i,j])                       
        exit()
        
        for j in range(1, ny):
            for i in range(1, nx):
                wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]

        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " wt[i,j]=",wt[i,j])    
                
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " wn[i,j]=",wn[i,j])                    
                
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                print("i,j=",i, " ", j, " r[i,j]=",r[i,j])                       
        exit()        
        """
        
        bc2(nx,ny,dx,dy,wn,sn)
        #print("sn[1,1]= ",sn[1,1] )
        
      

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)
        #cal_stream_function_DEBUG(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)
        #print("sn[1,1]= ",sn[1,1] )
        #exit()

        rms[k] = 0.0
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                ds = sn[i,j] - sp[i,j]
                #print(i," ", j, " ", ds)
                rms[k] = rms[k] + ds * ds

        rms[k] = np.sqrt( rms[k] / ( (nx+1)*(ny+1) ) )
        print(k, " ", rms[k])
        
def numericalTEST(nx,ny,nt,dx,dy,dt,re,wn,sn,rms,MaxIt,MaxErr,beta):
    wt = np.zeros((nx+1,ny+1)) # temporary array during RK3 integration
    r  = np.zeros((nx+1,ny+1)) # right hand side
    sp = np.zeros((nx+1,ny+1)) # old streamfunction
    
    for k in range(0, nt):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                sp[i,j] = sn[i,j]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,sn,r)
        
        for j in range(1, ny):
            for i in range(1, nx):
                wt[i,j] = wn[i,j] + dt*r[i,j]

        bc2(nx,ny,dx,dy,wt,sn)        

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for j in range(1, ny):
            for i in range(1, nx):
                wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]

        bc2(nx,ny,dx,dy,wt,sn)
        
        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)        
        
        for j in range(1, ny):
            for i in range(1, nx):
                wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]

        bc2(nx,ny,dx,dy,wn,sn)

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)
        #cal_stream_function_DEBUG(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)
        #print("sn[1,1]= ",sn[1,1] )
        #exit()

        rms[k] = 0.0
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                ds = sn[i,j] - sp[i,j]
                #print(i," ", j, " ", ds)
                rms[k] = rms[k] + ds * ds

        rms[k] = np.sqrt( rms[k] / ( (nx+1)*(ny+1) ) )
        print(k, " ", rms[k])        
        
def numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms,MaxIt,MaxErr,beta):
    wt = np.zeros((nx+1,ny+1)) # temporary array during RK3 integration
    r  = np.zeros((nx+1,ny+1)) # right hand side
    sp = np.zeros((nx+1,ny+1)) # old streamfunction
    
    for k in range(0, nt):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                sp[i,j] = sn[i,j]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,sn,r)
        
        for j in range(1, ny):
            for i in range(1, nx):
                wt[i,j] = wn[i,j] + dt*r[i,j]

        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for j in range(1, ny):
            for i in range(1, nx):
                wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]

        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r)
        
        for j in range(1, ny):
            for i in range(1, nx):
                wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]

        bc2(nx,ny,dx,dy,wn,sn)

        # compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta)

        rms[k] = 0.0
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                ds = sn[i,j] - sp[i,j]
                rms[k] = rms[k] + ds * ds

        rms[k] = np.sqrt( rms[k] / ( (nx+1)*(ny+1) ) )
        print(k, " ", rms[k])        

#nx = 64
#ny = 64

nx = 4
ny = 4

Visc=0.1;
dt = 0.001
#tf = 10.0
#tf = 1.0
#tf = 10.0
#tf = 0.001
tf = 0.01
nt = int(tf/dt)
re = 100.0

# parametersforSOR iteration
MaxIt=100;
Beta=1.5;
#MaxErr=0.001;
MaxErr=1.0e-6;

wn = np.zeros((nx+1,ny+1));
sn = np.zeros((nx+1,ny+1));

h=1.0/(nx);
t=0.0;

uwall = 1.0
x = np.zeros(nx+1)
y = np.zeros(ny+1)
rms = np.zeros(nt)

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

for j in range(1, ny+1):
    for i in range(0, nx+1):  
        wn[i,j] = 0.0 # initial condition
        sn[i,j] = 0.0 # initial streamfunction    

time = 0.0

#numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms,MaxIt,MaxErr,Beta)
numericalTEST(nx,ny,nt,dx,dy,dt,re,wn,sn,rms,MaxIt,MaxErr,Beta)

residual_plot = open("residual_plot.txt", "w");

for n in range(nt):
    residual_plot.write("{0} {1}\n" .format(n+1,rms[n]));
    
residual_plot.close()
  
# Creating 2-D grid of features
[X, Y] = np.meshgrid(x, y, indexing='ij')

fig = plt.figure("OneFLOW-CFD ", figsize=(10, 6), dpi=100)
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

cs = ax1.contourf(X, Y, wn, levels=60, cmap="jet")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_title("Vorticity field")

cs = ax2.contourf(X, Y, sn, levels=60, cmap="jet")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_title("Streamfunction")

plt.show()