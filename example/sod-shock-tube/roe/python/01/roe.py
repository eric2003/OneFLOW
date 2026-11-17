import numpy as np
import matplotlib.pyplot as plt

class MyPlot:
    def __init__( self ):
        self.x = []
        self.r = []
        self.u = []
        self.m = []
        self.p = []
    
    def AddData( self, xm, rm, um, mm, pm ):
        self.x.append( xm )
        self.r.append( rm )
        self.u.append( um )
        self.m.append( mm )
        self.p.append( pm )
        
    def ReadData( self ):
        with open('sod_theory.plt', 'r') as f:
            for index, line in enumerate(f):
                words = line.strip().split()
                self.x.append( float(words[0]) )
                self.r.append( float(words[1]) )
                self.u.append( float(words[2]) )
                self.m.append( float(words[3]) )
                self.p.append( float(words[4]) )
        self.ComputeEnergy()
                
    def ComputeEnergy( self ):
        num = len(self.x)
        self.e = np.zeros( num )
        print("self.e.len=", len(self.e))
        gama = 1.4
        for i in range(0, num ):
            um = self.u[i]
            rm = self.r[i]
            pm = self.p[i]
            self.e[i] = (1.0/(gama-1.0))* pm/rm + 0.5 * ( um * um )
    
    def PlotTheory( self ):
        plt.figure("Exact solution for the Sod's shock-tube problem", figsize=(10, 8), dpi=100)
        plt.subplot(2, 2, 1)
        plt.plot(self.x, self.r, linewidth=1.0, label="density")
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho$")
        plt.legend()
        
        plt.subplot(2, 2, 2)
        plt.plot(self.x, self.u, linewidth=1.0, label="velocity")
        plt.xlabel("$x$")
        plt.ylabel(r"$u$")
        plt.legend()
        
        plt.subplot(2, 2, 3)
        plt.plot(self.x, self.m, linewidth=1.0, label="mach number")
        plt.xlabel("$x$")
        plt.ylabel(r"$m$")
        plt.legend()
        
        plt.subplot(2, 2, 4)
        plt.plot(self.x, self.p, linewidth=1.0, label="pressure")
        plt.xlabel("$x$")
        plt.ylabel(r"$p$")
        plt.legend()
        
        plt.tight_layout()
        
        plt.show()
        
    def Plot( self, x, q ):
        plt.figure("OneFLOW-CFD Solver", figsize=(6, 4), dpi=100)
        plt.subplot(2, 2, 1)
        plt.plot(x, q[:,0, 0], linewidth=1.0, label="t=")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$\rho$")
        plt.title(r"rho")
        plt.legend(loc='upper right', fontsize='6')
        
        plt.subplot(2, 2, 2)
        plt.plot(x, q[:,1, 0], linewidth=1.0, label="t=")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$u$")
        plt.title(r"u")
        plt.legend(loc='upper right', fontsize='6')
        
        plt.subplot(2, 2, 3)
        plt.plot(x, q[:,2, 0], linewidth=1.0, label="t=")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$e$")
        plt.title(r"e")
        plt.legend(loc='upper right', fontsize='6') 
        
        plt.tight_layout()
        plt.show()  

    def PlotCompare( self, x, q ):
        numPoints = len( q[:,0, 0] )
        print("numPoints=",numPoints)
        
        rr = np.zeros( numPoints )
        uu = np.zeros( numPoints )
        pp = np.zeros( numPoints )
        ee = np.zeros( numPoints )
        
        gama = 1.4
        for i in range( 0, numPoints ):
            rho  = q[ i, 0, 0 ]
            rhou = q[ i, 1, 0 ]
            rhoe = q[ i, 2, 0 ]
            rr[i] = rho
            uu[i] = rhou / rho
            ee[i] = rhoe / rho
            pp[i] = ( gama - 1.0 ) * rhoe - 0.5 * rho * ( uu[i] * uu[i] ) 
          
        
        plt.figure("Exact solution for the Sod's shock-tube problem", figsize=(10, 8), dpi=100)
        plt.subplot(2, 2, 1)
        plt.plot(self.x, self.r, linewidth=1.0, label="theory")
        #plt.plot(x, q[:,0, 0], linewidth=1.0, label="OneFLOW-CFD")
        plt.plot(x, rr, linewidth=1.0, label="OneFLOW-CFD")
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho$")
        plt.legend()
        
        plt.subplot(2, 2, 2)
        plt.plot(self.x, self.u, linewidth=1.0, label="velocity")
        #plt.plot(x, q[:,1, 0], linewidth=1.0, label="OneFLOW-CFD")
        plt.plot(x, uu, linewidth=1.0, label="OneFLOW-CFD")
        plt.xlabel("$x$")
        plt.ylabel(r"$u$")
        plt.legend()
        
        plt.subplot(2, 2, 3)
        plt.plot(self.x, self.e, linewidth=1.0, label="theory")
        plt.plot(x, ee, linewidth=1.0, label="OneFLOW-CFD")
        plt.xlabel("$x$")
        plt.ylabel(r"$E$")
        plt.legend()
        
        plt.subplot(2, 2, 4)
        plt.plot(self.x, self.p, linewidth=1.0, label="pressure")
        #plt.plot(x, q[:,1, 0], linewidth=1.0, label="OneFLOW-CFD")
        plt.plot(x, pp, linewidth=1.0, label="OneFLOW-CFD")
        plt.xlabel("$x$")
        plt.ylabel(r"$p$")
        plt.legend()
        
        plt.tight_layout()
        
        plt.show()        

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
def fluxes(nx, gamma, q, f ):
    #i=0,1,...,nx
    for i in range(0, nx + 1):
        p = ( gamma - 1.0 ) * ( q[i,2] - 0.5 * q[i,1] * q[i,1] / q[i,0] )
        f[i,0] = q[i,1]
        f[i,1] = q[i,1]*q[i,1]/q[i,0] + p
        f[i,2] = q[i,1]*q[i,2]/q[i,0] + p * q[i,1]/q[i,0]
    return
    
    
#-----------------------------------------------------------------------------#
# Riemann solver: Roe's approximate Riemann solver
#-----------------------------------------------------------------------------#
def roe( nx, gamma, qL, qR, fL, fR, f ):
    dd = np.zeros( 3 )
    dF = np.zeros( 3 )
    dQ = np.zeros( 3 )
    gm = gamma - 1.0
    
    #i=0,1,...,nx
    for i in range(0, nx+1):
        #Left and right states:
        rhoL = qL[i,0]
        uL   = qL[i,1]/rhoL
        eL   = qL[i,2]/rhoL
        pL   = gm * ( rhoL * eL - 0.5 * rhoL * ( uL * uL ) )
        hL   = eL + pL / rhoL
    
        rhoR = qR[i,0]
        uR   = qR[i,1]/rhoR
        eR   = qR[i,2]/rhoR
        pR   = gm * ( rhoR * eR - 0.5 * rhoR * ( uR * uR ) )
        hR   = eR + pR / rhoR
    
        alpha = 1.0/(np.sqrt(abs(rhoL)) + np.sqrt(abs(rhoR)))
    
        uu = ( np.sqrt(abs(rhoL)) * uL + np.sqrt(abs(rhoR)) * uR ) * alpha
        hh = ( np.sqrt(abs(rhoL)) * hL + np.sqrt(abs(rhoR)) * hR ) * alpha
        aa = np.sqrt( abs( gm * ( hh - 0.5*uu*uu ) ) )
    
        D11 = abs( uu )
        D22 = abs( uu + aa )
        D33 = abs( uu - aa )
    
        beta = 0.5/( aa * aa )
        phi2 = 0.5 * gm * uu * uu
    
        #Right eigenvector matrix
        R11, R21, R31 = 1.0, uu, phi2/gm
        R12, R22, R32 = beta, beta*(uu + aa), beta*(hh + uu*aa)
        R13, R23, R33 = beta, beta*(uu - aa), beta*(hh - uu*aa)
    
        #Left eigenvector matrix
        L11, L12, L13 = 1.0-phi2/(aa*aa), gm*uu/(aa*aa), -gm/(aa*aa)
        L21, L22, L23 = phi2 - uu*aa, - gm*uu + aa, gm
        L31, L32, L33 = phi2 + uu*aa, - gm*uu - aa, gm
    
        for m in range(0, 3):
            dQ[m] = qR[i,m] - qL[i,m]
    
        dd[0] = D11 * ( L11 * dQ[0] + L12 * dQ[1] + L13 * dQ[2] )
        dd[1] = D22 * ( L21 * dQ[0] + L22 * dQ[1] + L23 * dQ[2] )
        dd[2] = D33 * ( L31 * dQ[0] + L32 * dQ[1] + L33 * dQ[2] )
    
        dF[0] = R11 * dd[0] + R12 * dd[1] + R13 * dd[2]
        dF[1] = R21 * dd[0] + R22 * dd[1] + R23 * dd[2]
        dF[2] = R31 * dd[0] + R32 * dd[1] + R33 * dd[2]
    
        for m in range(0, 3):
            f[i,m] = 0.5 * ( fR[i,m] + fL[i,m] ) - 0.5 * dF[m]

#-----------------------------------------------------------------------------#
# Calculate right hand side terms of the Euler equations
#-----------------------------------------------------------------------------#
def rhs( nx, dx, gamma, q, r ):
    qL  = np.zeros( ( nx+1, 3 ) )
    qR  = np.zeros( ( nx+1, 3 ) )
    
    # left and right side fluxes at the interface
    fL = np.zeros( ( nx+1, 3 ) )
    fR = np.zeros( ( nx+1, 3 ) )
    
    #fluxes at the interface    
    f = np.zeros( ( nx+1, 3 ) )
    
    # WENO Reconstruction
    wenoL( nx, q, qL )
    wenoR( nx, q, qR )
    
    # Computing fluxes
    fluxes( nx, gamma, qL, fL )
    fluxes( nx, gamma, qR, fR )    

    # compute Riemann solver using Roe scheme(flux at interface)
    roe( nx, gamma, qL, qR, fL, fR, f )

    # RHS
    for i in range(0, nx):
        for m in range(0, 3):
            r[i,m] = -( f[i+1,m] - f[i,m] ) / dx


#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 5th-order Compact WENO scheme for spatial terms
#-----------------------------------------------------------------------------#
def numerical( nx, ns, nt, dx, dt, u ):
    x  = np.zeros( nx )
    qn = np.zeros( ( nx, 3 ) ) # numerical solsution at every time step
    qt = np.zeros( ( nx, 3 ) ) # temporary array during RK3 integration
    r  = np.zeros( ( nx, 3 ) )
    
    k = 0 # record index
    freq = int( nt / ns )
    print("freq=",freq)
    
    gamma = 1.4 # specific gas ratio    
    
    # Sod's Riemann problem
    # Left side
    rhoL = 1.0
    uL   = 0.0
    pL   = 1.0
    # Right side
    rhoR = 0.125
    uR   = 0.0
    pR   = 0.1

    # nodal storage location (grid)
    #i=0,1,...,nx-1
    for i in range(0, nx):
        x[i] = dx * ( i + 0.5 )
        
    xc = 0.5 # seperator location
    #i=0,1,...,nx-1
    for i in range(0, nx):
        if ( x[i] > xc ):
            rho = rhoR
            u = uR
            p = pR
        else:
            rho = rhoL
            u = uL
            p = pL
    
        e = p / ( rho * ( gamma - 1.0 ) ) + 0.5 * u * u
    
        #conservative variables
        qn[i,0] = rho
        qn[i,1] = rho * u
        qn[i,2] = rho * e
        
    #i=0,1,...,nx-1
    for i in range(0, nx):
        #m=0,1,2
        for m in range(0, 3):
            q[i,m,k] = qn[i,m] # store solution at t=0
            
    mplot = MyPlot()
    mplot.ReadData()
    mplot.PlotTheory()
    mplot.Plot(x, q)
    mplot.PlotCompare(x, q)
            
    exit(0)                    
    # TVD RK3 for time integration
    #n=1,2,...,nt
    for n in range(1, nt+1): # time step
        rhs( nx, dx, gamma, qn, r)

        for i in range(0, nx):
            for m in range(0, 3):
                qt[i,m] = qn[i,m] + dt*r[i,m]
                
        rhs( nx, dx, gamma, qt, r)   

        for i in range(0, nx):
            for m in range(0, 3):
                qt[i,m] = 0.75*qn[i,m] + 0.25*qt[i,m] + 0.25*dt*r[i,m]
            
        rhs( nx, dx, gamma, qt, r )

        for i in range(0, nx):
            for m in range(0, 3):
                qn[i,m] = (1.0/3.0)*qn[i,m] + (2.0/3.0)*qt[i,m] + (2.0/3.0)*dt*r[i,m]
            
        if ( n% freq == 0 ):
            k = k+1  
            for i in range(0, nx):
                for m in range(0, 3):
                   q[i,m,k] = qn[i,m]
            print("k=",k,"ns=",ns)

nx = 256
ns = 20
dt = 0.0001
tm = 0.20

dx = 1.0 / nx
nt = int( tm / dt )
ds = tm / ns

print("nx=",nx)
print("ns=",ns)
print("dt=",dt)
print("tm=",tm)
print("dx=",dx)
print("nt=",nt)

q  = np.zeros( ( nx, 3, ns+1 ) )
numerical(nx,ns,nt,dx,dt,q)

x = np.linspace(0,1, num=nx+1) 

plt.figure("OneFLOW-CFD Solver", figsize=(6, 4), dpi=100)
plt.plot(x, q[:,0, 0], linewidth=1.0, label="t=")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Inviscid Burgers Equation: WENO5+Rusanov Scheme+Periodic BC")
plt.legend(loc='upper right', fontsize='6')
plt.show()    
