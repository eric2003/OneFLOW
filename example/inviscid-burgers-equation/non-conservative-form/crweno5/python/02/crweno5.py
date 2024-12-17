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
    
#-----------------------------------------------------------------------------#
# Solution to tridigonal system using Thomas algorithm
#-----------------------------------------------------------------------------#
def thomas_algorithm(a, b, c, d, n):
    c_prime = [0] * n
    d_prime = [0] * n
    x = [0] * n
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n):
        coef = 1.0 / ( b[i] - a[i] * c_prime[i-1] )
        c_prime[i] = c[i] * coef
        d_prime[i] = ( d[i] - a[i] * d_prime[i-1] ) * coef

    x[n-1] = d_prime[n-1]
    

    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    
    return x     
    
#---------------------------------------------------------------------------#
#nonlinear weights for upwind direction
#---------------------------------------------------------------------------#
def crwcL(v1,v2,v3,v4,v5):
    eps = 1.0e-6

    # smoothness indicators
    s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

    # computing nonlinear weights w1,w2,w3
    c1 = 2.0e-1 / ( (eps+s1)**2 )
    c2 = 5.0e-1 / ( (eps+s2)**2 )
    c3 = 3.0e-1 / ( (eps+s3)**2 )

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)
    
    a1 = (2.0*w1 + w2)/3.0
    a2 = (w1 + 2.0*w2 + 2.0*w3)/3.0
    a3 = w3/3.0

    b1 = w1/6.0
    b2 = (5.0*w1 + 5.0*w2 + w3)/6.0
    b3 = (w2 + 5.0*w3)/6.0

    return a1,a2,a3,b1,b2,b3    


#---------------------------------------------------------------------------#
#nonlinear weights for downwind direction
#---------------------------------------------------------------------------#
def crwcR(v1,v2,v3,v4,v5):
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

    c1 = 3.0e-1/(eps+s1)**2
    c2 = 5.0e-1/(eps+s2)**2
    c3 = 2.0e-1/(eps+s3)**2

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    a1 = w1/3.0
    a2 = (w3 + 2.0*w2 + 2.0*w1)/3.0
    a3 = (2.0*w3 + w2)/3.0

    b1 = (w2 + 5.0*w1)/6.0
    b2 = (5.0*w3 + 5.0*w2 + w1)/6.0
    b3 = w3/6.0

    return a1,a2,a3,b1,b2,b3
    
#-----------------------------------------------------------------------------#
# CRWENO reconstruction ofr upwind direction (positive and left to right)
# u(i): solution values at finite difference grid nodes i = 0,1,...,N
# f(j): reconstructed values at nodes j = i+1/2; j = 0,1,...,N-1
#-----------------------------------------------------------------------------#
def crwenoL(nx,u,f):
    a = np.zeros(nx)
    b = np.zeros(nx)
    c = np.zeros(nx)
    r = np.zeros(nx)

    i = 0    
    b[i] = 2.0/3.0
    c[i] = 1.0/3.0
    r[i] = (u[i] + 5.0*u[i+1])/6.0
    
    i = 1
    v1 = 2.0*u[i-1] - u[i]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = u[i+2]
    
    a1,a2,a3,b1,b2,b3 = crwcL(v1,v2,v3,v4,v5)
    a[i] = a1
    b[i] = a2
    c[i] = a3
    r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]
    
    for i in range(2, nx-1):
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i  ]
        v4 = u[i+1]
        v5 = u[i+2]
        a1,a2,a3,b1,b2,b3 = crwcL(v1,v2,v3,v4,v5)
        a[i] = a1
        b[i] = a2
        c[i] = a3
        r[i] = b1*u[i-1] + b2*u[i] + b3*u[i+1]
        
    #i = nx-1
    #a[i] = 1.0/3.0
    #b[i] = 2.0/3.0
    #r[i] = (5.0*u[i] + u[i+1])/6.0        

    i = nx-1
    a[i] = 2.0/3.0
    b[i] = 1.0/3.0
    r[i] = (u[i] + 5.0*u[i+1])/6.0
    
    f[:] = thomas_algorithm( a, b, c, r, nx )
    

#-----------------------------------------------------------------------------#
# CRWENO reconstruction for downwind direction (negative; right to left)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
#-----------------------------------------------------------------------------#
def crwenoR(nx,u,f):
    a = np.zeros(nx)
    b = np.zeros(nx)
    c = np.zeros(nx)
    r = np.zeros(nx)

    #i = 1
    #b[i-1] = 2.0/3.0
    #c[i-1] = 1.0/3.0
    #r[i-1] = (u[i-1] + 5.0*u[i])/6.0    

    i = 1
    b[i-1] = 1.0/3.0
    c[i-1] = 2.0/3.0
    r[i-1] = (5.0*u[i-1] + u[i])/6.0    

    for i in range(2, nx-1):
        v1 = u[i-2]
        v2 = u[i-1]
        v3 = u[i]
        v4 = u[i+1]
        v5 = u[i+2]

        a1,a2,a3,b1,b2,b3 = crwcR(v1,v2,v3,v4,v5)
        a[i-1] = a1
        b[i-1] = a2
        c[i-1] = a3
        r[i-1] = b1*u[i-1] + b2*u[i] + b3*u[i+1]
        
    i = nx-1
    v1 = u[i-2]
    v2 = u[i-1]
    v3 = u[i]
    v4 = u[i+1]
    v5 = 2.0*u[i+1] - u[i]    

    a1,a2,a3,b1,b2,b3 = crwcR(v1,v2,v3,v4,v5)
    a[i-1] = a1
    b[i-1] = a2
    c[i-1] = a3
    r[i-1] = b1*u[i-1] + b2*u[i] + b3*u[i+1]        

    i = nx
    a[i-1] = 1.0/3.0
    b[i-1] = 2.0/3.0
    r[i-1] = (5.0*u[i-1] + u[i])/6.0

    f[:] = thomas_algorithm( a, b, c, r, nx )


#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -u∂u/∂x
#-----------------------------------------------------------------------------#
def rhs( nx, dx, u, r ):
    uL = np.zeros( nx  )
    uR = np.zeros( nx  )

    crwenoL( nx, u, uL )
    crwenoR( nx, u, uR )

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
            
        rhs( nx, dx, ut, r )

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
plt.title("Inviscid Burgers Equation: CRWENO-5 Scheme+Dirichlet BC")
plt.legend(loc='upper right', fontsize='6')
plt.show()    
