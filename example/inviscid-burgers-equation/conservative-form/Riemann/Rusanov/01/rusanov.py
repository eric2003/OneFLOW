import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
def compute_l2norm(nx,r):
    rms = 0.0
    for i in range(1, nx):
        rms += r[i] * r[i]
    rms = np.sqrt( rms / ( ( nx - 1 ) ) )
    return rms
    
#---------------------------------------------------------------------------#
#nonlinear weights for upwind direction
#---------------------------------------------------------------------------#
def wcL(v1,v2,v3,v4,v5):
    eps = 1.0e-6

    # smoothness indicators
    s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

    # computing nonlinear weights w1,w2,w3
    c1 = 1.0e-1 / ( (eps+s1)**2 )
    c2 = 6.0e-1 / ( (eps+s2)**2 )
    c3 = 3.0e-1 / ( (eps+s3)**2 )

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    # candiate stencils
    q1 = v1/3.0 - 7.0/6.0*v2 + 11.0/6.0*v3
    q2 =-v2/6.0 + 5.0/6.0*v3 + v4/3.0
    q3 = v3/3.0 + 5.0/6.0*v4 - v5/6.0

    # reconstructed value at interface
    f = (w1*q1 + w2*q2 + w3*q3)

    return f


#---------------------------------------------------------------------------#
#nonlinear weights for downwind direction
#---------------------------------------------------------------------------#
def wcR(v1,v2,v3,v4,v5):
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

    c1 = 3.0e-1/(eps+s1)**2
    c2 = 6.0e-1/(eps+s2)**2
    c3 = 1.0e-1/(eps+s3)**2

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    # candiate stencils
    q1 =-v1/6.0      + 5.0/6.0*v2 + v3/3.0
    q2 = v2/3.0      + 5.0/6.0*v3 - v4/6.0
    q3 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0

    # reconstructed value at interface
    f = (w1*q1 + w2*q2 + w3*q3)

    return f
    
#-----------------------------------------------------------------------------#
# WENO reconstruction for upwind direction (positive; left to right)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i+1/2; j = 1,...,N
#-----------------------------------------------------------------------------#
def wenoL(nx,u,f):
    i = -1
    v1 = u[nx-3]
    v2 = u[nx-2]
    v3 = u[nx-1]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i+1] = wcL(v1,v2,v3,v4,v5)

    i = 0
    v1 = u[nx-2]
    v2 = u[nx-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i+1] = wcL(v1,v2,v3,v4,v5)
    
    i = 1
    v1 = u[nx-1]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i+1] = wcL(v1,v2,v3,v4,v5)
    
    #i=2,3,...,nx-4,nx-3
    for i in range(2, nx-2):
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i]
        v4 = u[i+1]
        v5 = u[i+2]
        f[i+1] = wcL(v1,v2,v3,v4,v5)

    i = nx-2
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i  ]
    v4 = u[i+1]
    v5 = u[0]
    f[i+1] = wcL(v1,v2,v3,v4,v5) 

    i = nx-1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[0]
    v5 = u[1]
    f[i+1] = wcL(v1,v2,v3,v4,v5) 
    
#-----------------------------------------------------------------------------#
# WENO reconstruction for downwind direction (negative; right to left)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
#-----------------------------------------------------------------------------#
def wenoR(nx,u,f):
    i = 0
    v1 = u[nx-2]
    v2 = u[nx-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i] = wcR(v1,v2,v3,v4,v5)
    
    i = 1
    v1 = u[nx-1]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i] = wcR(v1,v2,v3,v4,v5)
    
    #i=2,3,...,nx-4,nx-3
    for i in range(2, nx-2):
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i  ]
        v4 = u[i+1]
        v5 = u[i+2]
        f[i] = wcR(v1,v2,v3,v4,v5)
        
    i = nx-2
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[0]
    f[i] = wcR(v1,v2,v3,v4,v5)

    i = nx-1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[0]
    v5 = u[1]
    f[i] = wcR(v1,v2,v3,v4,v5)

    i = nx
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[0]
    v4 = u[1]
    v5 = u[2]
    f[i] = wcR(v1,v2,v3,v4,v5)  
    
#-----------------------------------------------------------------------------#
# Calculate fluxes
#-----------------------------------------------------------------------------#
def burgers_fluxes( nx, u, f ):
    #i=0,1,...,nx
    for i in range(0, nx + 1):
        f[i] = 0.5 * u[i] * u[i]
    return
    
#-----------------------------------------------------------------------------#
# Riemann solver: Rusanov
#-----------------------------------------------------------------------------#
def rusanov( nx, u, uL, uR, fL, fR, f ):

    # propagation speed
    ps = np.zeros( nx + 1  )
    
    ps[0] = max( abs(u[0]), abs(u[nx-1]) )
    #i=1,2,...,nx-1
    for i in range(1, nx):
        ps[i] = max( abs(u[i]), abs(u[i-1]) )

    ps[nx] = max(abs(u[0]), abs(u[nx-1]))

    # Interface fluxes (Rusanov)
    for i in range(0, nx+1):
        f[i] = 0.5 * ( fR[i] + fL[i] ) - 0.5 * ps[i] * ( uR[i] - uL[i] )
    return


#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -u∂u/∂x
#-----------------------------------------------------------------------------#
def rhs( nx, dx, u, r ):
    uL = np.zeros( nx + 1 )
    uR = np.zeros( nx + 1 )

    # left and right side fluxes at the interface
    fL = np.zeros( nx + 1 )
    fR = np.zeros( nx + 1 )
    
    #fluxes at the interface    
    f = np.zeros( nx + 1 )
    
    # WENO Reconstruction
    wenoL( nx, u, uL )
    wenoR( nx, u, uR )
    
    # Computing fluxes
    burgers_fluxes( nx, uL, fL )
    burgers_fluxes( nx, uR, fR )

    # compute Riemann solver (flux at interface)
    rusanov( nx, u, uL, uR, fL, fR, f )

    # RHS
    for i in range(0, nx):
        r[i] = - ( f[i+1] - f[i] ) / dx


#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 5th-order Compact WENO scheme for spatial terms
#-----------------------------------------------------------------------------#
def numerical( nx, ns, nt, dx, dt, u ):
    x  = np.zeros(nx+1)
    un = np.zeros(nx+1) # numerical solsution at every time step
    ut = np.zeros(nx+1) # temporary array during RK3 integration
    r  = np.zeros(nx)

    k = 0 # record index
    freq = int( nt / ns )
    print("freq=",freq)

    for i in range(0, nx):
        x[i] = dx * ( i + 0.5 )
        un[i] = np.sin( 2.0 * np.pi * x[i] )
        u[i,k] = un[i] # store solution at t=0
        

    for j in range(1, nt+1):
        rhs( nx, dx, un, r )

        for i in range(0, nx):
            ut[i] = un[i] + dt*r[i]
            
           
        rhs( nx, dx, ut, r )

        for i in range(0, nx):
            ut[i] = 0.75*un[i] + 0.25*ut[i] + 0.25*dt*r[i]
            
           
        rhs( nx, dx, ut, r )

        for i in range(0, nx):
            un[i] = (1.0/3.0)*un[i] + (2.0/3.0)*ut[i] + (2.0/3.0)*dt*r[i]
            
            
        if ( j% freq == 0 ):
            k = k+1        
            u[:,k] = un[:]
            print("k=",k,"ns=",ns)

nx = 200
ns = 10
dt = 0.0001
tm = 0.25

dx = 1.0 / nx
nt = int( tm / dt )
ds = tm / ns

print("nx=",nx)
print("ns=",ns)
print("dt=",dt)
print("tm=",tm)
print("dx=",dx)
print("nt=",nt)

u = np.zeros( (nx+1, ns+1 ) )
numerical( nx, ns, nt, dx, dt, u )

x = np.linspace(0,1, num=nx+1) 

plt.figure("OneFLOW-CFD Solver", figsize=(6, 4), dpi=100)
for k in range(0, ns+1):
    plt.plot(x, u[:,k], linewidth=1.0, label="t="+format(tm*k/ns, ".4f"))
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Inviscid Burgers Equation: WENO5+Rusanov Scheme+Periodic BC")
plt.legend(loc='upper right', fontsize='6')
plt.show()    
