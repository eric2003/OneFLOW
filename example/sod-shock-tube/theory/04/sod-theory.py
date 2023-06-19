import numpy as np
import matplotlib.pyplot as plt

class MyPlot:
    def __init__( self):
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
    
    def Plot( self ):
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

def sp2p1( gam, p1, a1, p4, a4, tol ):
    #Uses Newton-secant method to iterate on eqn 7.94 (Anderson, 
    #1984) to fine p2p1 across moving shock wave.

    gm1 = gam - 1.0;
    gp1 = gam + 1.0;

    #Initialize p2p1 for starting guess
    p2p1m = 0.9 * p4 / p1;
    t1 = - 2.0 * gam / gm1;

    t2 = gm1 * ( a1 / a4 ) * ( p2p1m - 1.0 );
    t3 = 2.0 * gam * ( 2.0 * gam + gp1 * ( p2p1m - 1.0 ) );
    fm = p4 / p1 - p2p1m * pow( 1.0 - t2 / np.sqrt( t3 ), t1 );

    #Perturb p2p1
    p2p1 = 0.95 * p2p1m;

    #Begin iteration
    iter = 0;
    itmax = 20;

    while True:
        iter = iter + 1

        t2 = gm1 * ( a1 / a4 ) * ( p2p1 - 1.0 )
        t3 = 2.0 * gam * ( 2.0 * gam + gp1 * ( p2p1 - 1.0 ) )

        f = p4 / p1 - p2p1 * pow( 1.0 - t2 / np.sqrt( t3 ), t1 );

        print( "iter, p2p1, f: ", iter, " ", p2p1, " ", f )

        if abs( f ) <= tol or iter >= itmax:
            break

        p2p1n = p2p1 - f * ( p2p1 - p2p1m ) / ( f - fm );
        p2p1m = p2p1;
        fm = f;
        p2p1 = p2p1n;

    #Check to see if maximum iterations reached
    iterr = 0;
    if iter < itmax:
        iterr = 1
    return p2p1, iterr
    
def Sod_Theory():
    #  Calcs the exact solution for the Sod's shock-tube problem.
    #  See J.D. Anderson, Modern Compressible Flow (1984) for details.
    #
    #     --------------------------------------------
    #     |                     |                    |
    #     |    p4, r4, u4       |     p1, r1, u1     |   tm = 0.0
    #     |                     |                    |
    #     --------------------------------------------
    #     xl                    xd                   xr
    #                p4>p1
    #
    #
    #     --------------------------------------------
    #     |      |      |           | up       | W   |
    #     |   4  |<-----|    3    --|->  2   --|-->  |
    #     |      |      |           |          |     |
    #     --------------------------------------------
    #     xl    expansion          slip     shock    xr
    #
    #-----------------------------------------------------------------

    #Set constants.
    gam = 1.4;
    gm1 = gam - 1.0;
    gp1 = gam + 1.0;

    #Set initial states (non-dimensional).
    p4 = 1.0;
    r4 = 1.0;
    u4 = 0.0;

    p1 = 0.1;
    r1 = 0.125;
    u1 = 0.0;

    tol = 1.0E-05;
    tm = 0.20;

    #Set dimensions of shocktube.
    xl = 0.0;
    xr = 1.0;
    xd = 0.5;

    #Calc acoustic velocities.
    a1 = np.sqrt( gam * p1 / r1 );
    a4 = np.sqrt( gam * p4 / r4 );

    #Use a Newton-secant iteration to compute p2p1.
    p2p1, iterr = sp2p1( gam, p1, a1, p4, a4, tol );

    print( "p2p1 = ", p2p1 )
    
    t2t1 = p2p1 * ( gp1 / gm1 + p2p1 ) / ( 1.0 + gp1 * p2p1 / gm1 );

    r2r1 = ( 1.0 + gp1 * p2p1 / gm1 ) / ( gp1 / gm1 + p2p1 );

    #shock-wave speed.
    wsp = a1 * np.sqrt( gp1 * ( p2p1 - 1.0 ) / ( 2.0 * gam ) + 1.0 );

    #Shock location.
    xs = xd + wsp * tm;

    #State 2.
    p2 = p2p1 * p1;
    r2 = r2r1 * r1;

    #State 3.
    p3 = p2;

    #Isentropic between 3 and 4.
    r3 = r4 * pow( p3 / p4, 1.0 / gam );
    a3 = np.sqrt( gam * p3 / r3 );

    #Speed of contact discontinuity.
    up = 2.0 * a4 * ( 1.0 - pow( p2 / p4, 0.5 * gm1 / gam ) )/ gm1;
    u2 = up;
    u3 = up;

    #Mach numbers.
    rmach1 = u1 / np.sqrt( gam * p1 / r1 );
    rmach2 = u2 / np.sqrt( gam * p2 / r2 );
    rmach3 = u3 / np.sqrt( gam * p3 / r3 );
    rmach4 = u4 / np.sqrt( gam * p4 / r4 );

    #Location of contact discontinuity.
    xc = xd + up * tm;

    #Location of expansion region.
    xhead = xd + ( u4 - a4 ) * tm;
    xtail = xd + ( u3 - a3 ) * tm;

    #Write out some data.
    print()
    print( "gamma             = ", gam)
    print( "diaphram location = ", xd )
    print( "time              = ", tm )
    print()
    print( "(1) p1 = ", p1 )
    print( "    r1 = ", r1 )
    print( "    u1 = ", u1 )
    print( "    m1 = ", rmach1 )
    print()
    print( "Shock speed    = ", wsp )
    print( "Shock location = ", xs )
    print()
    print( "(2) p2 = ", p2 )
    print( "    r2 = ", r2 )
    print( "    u2 = ", u2 )
    print( "    m2 = ", rmach2 )
    print()
    print( "Contact discontinuity speed    = ", up )
    print( "Contact discontinuity location = ", xc )
    print()
    print( "(3) p3 = ", p3 )
    print( "    r3 = ", r3 )
    print( "    u3 = ", u3 )
    print( "    m3 = ", rmach3 )
    print()
    print( "Expansion region head = ", xhead )
    print( "Expansion region tail = ", xtail )
    print()
    print( "(4) p4 = ", p4  )
    print( "    r4 = ", r4  )
    print( "    u4 = ", u4  )
    print( "    m4 = ", rmach4 )
    
    #Write out to files.
    output = open("sod_theory.plt", "w");
    nxp = 21;
    nNodes = nxp + 8;
    
    mplot = MyPlot();
    
    mplot.AddData(xl, r4, u4, rmach4, p4)
    mplot.AddData(xhead, r4, u4, rmach4, p4)

    for n in range(1, nxp+1):
        xx = xhead + ( xtail - xhead )  * n / ( nxp + 1.0 );
        ux = u4 + ( u3 - u4 ) * ( xx - xhead ) / ( xtail - xhead );
        px = p4 * pow( 1.0 - 0.5 * gm1 * ( ux / a4 ), 2.0 * gam / gm1 );
        rx = r4 * pow( 1.0 - 0.5 * gm1 * ( ux / a4 ), 2.0 / gm1 );
        mx = ux / np.sqrt( gam * px / rx );
        mplot.AddData(xx, rx, ux, mx, px)
        
    mplot.AddData(xtail, r3, u3, rmach3, p3)
    mplot.AddData(xc   , r3, u3, rmach3, p3)
    mplot.AddData(xc   , r2, u2, rmach2, p2)
    mplot.AddData(xs   , r2, u2, rmach2, p2)
    mplot.AddData(xs   , r1, u1, rmach1, p1)
    mplot.AddData(xr   , r1, u1, rmach1, p1)
    
    for i in range(0, len(mplot.x) ):
        output.write( "{0}  {1}  {2}  {3}  {4}\n".format(mplot.x[i], mplot.r[i], mplot.u[i], mplot.m[i], mplot.p[i]) )
    output.close()
    
    mplot.Plot()

Sod_Theory()
