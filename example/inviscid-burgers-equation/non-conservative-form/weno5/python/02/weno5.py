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
    a = np.zeros(nx)
    b = np.zeros(nx)
    c = np.zeros(nx)
    r = np.zeros(nx)

    i = 0
    v1 = 3.0*u[i] - 2.0*u[i+1]
    v2 = 2.0*u[i] - u[i+1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i] = wcL(v1,v2,v3,v4,v5)

    i = 1
    v1 = 2.0*u[i-1] - u[i]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i] = wcL(v1,v2,v3,v4,v5)

    for i in range(2, nx-1):
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i]
        v4 = u[i+1]
        v5 = u[i+2]
        f[i] = wcL(v1,v2,v3,v4,v5)

    i = nx-1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = 2.0*u[i+1]-u[i]
    f[i] = wcL(v1,v2,v3,v4,v5) 
    
#-----------------------------------------------------------------------------#
# CRWENO reconstruction for downwind direction (negative; right to left)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
#-----------------------------------------------------------------------------#
def wenoR(nx,u,f):
    i = 1
    v1 = 2.0*u[i-1] - u[i]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    f[i-1] = wcR(v1,v2,v3,v4,v5)
    
    for i in range(2, nx-1):
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i  ]
        v4 = u[i+1]
        v5 = u[i+2]
        f[i-1] = wcR(v1,v2,v3,v4,v5)

    i = nx-1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = 2.0*u[i+1] - u[i]
    f[i-1] = wcR(v1,v2,v3,v4,v5)

    i = nx
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = 2.0*u[i] - u[i-1]
    v5 = 3.0*u[i] - 2.0*u[i-1]
    f[i-1] = wcR(v1,v2,v3,v4,v5)    

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -u∂u/∂x
#-----------------------------------------------------------------------------#
def rhs( nx, dx, u, r ):
    uL = np.zeros( nx  )
    uR = np.zeros( nx  )

    wenoL( nx, u, uL )
    wenoR( nx, u, uR ) 

    for i in range(1, nx):
        if ( u[i] >= 0.0 ):
            r[i] = - u[i] * ( uL[i] - uL[i-1] ) / dx
        else:
            r[i] = - u[i] * ( uR[i] - uR[i-1] ) / dx

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

    for i in range(0, nx+1):
        x[i] = dx * ( i )
        un[i] = np.sin( 2.0 * np.pi * x[i] )
        u[i,k] = un[i] # store solution at t=0

    # dirichlet boundary condition
    u[0 ,k] = 0.0
    u[nx,k] = 0.0
    
    un[0 ] = 0.0
    un[nx] = 0.0

    # dirichlet boundary condition for temporary array
    ut[0 ] = 0.0
    ut[nx] = 0.0
    
    for j in range(1, nt+1):
        rhs( nx, dx, un, r )

        for i in range(1, nx):
            ut[i] = un[i] + dt*r[i]

        rhs( nx, dx, ut, r)

        for i in range(1, nx):
            ut[i] = 0.75*un[i] + 0.25*ut[i] + 0.25*dt*r[i]

        rhs( nx, dx, ut, r )

        for i in range(1, nx):
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
plt.title("Inviscid Burgers Equation: Non-Conservative Form-WENO-5 Scheme")
plt.legend(loc='upper right', fontsize='6')
plt.show()    


