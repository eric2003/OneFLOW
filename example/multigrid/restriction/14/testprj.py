import matplotlib.pyplot as plt
import numpy as np

class MG:
    def __init__(self):
        self.level = 0
        
    def allocate( self, nx, ny, dx, dy, level ):
        self.level = int(level)
        self.nx = int(nx)
        self.ny = int(ny)
        self.dx = dx
        self.dy = dy
        self.vh = np.zeros( ( nx + 1, ny + 1 ), dtype=np.float64)
        self.fh = np.zeros( ( nx + 1, ny + 1 ), dtype=np.float64)

class Mgdata:
    def __init__(self):
        self.nlevel = 1
        self.v1 = 2 # relaxation
        self.v2 = 2 # prolongation
        self.v3 = 2 # coarsest level           
    
    def allocate( self, nxm, nym, dxm, dym, nlevel ):
        self.nlevel = nlevel
        self.nxm = nxm
        self.nym = nym
        self.dxm = dxm
        self.dym = dym
        #allocate memory for grid size at different levels
        self.nx = np.zeros(nlevel, dtype=np.int32 )
        self.ny = np.zeros(nlevel, dtype=np.int32 )
        self.dx = np.zeros(nlevel, dtype=np.float64)
        self.dy = np.zeros(nlevel, dtype=np.float64)
        
        self.nx[0] = nxm
        self.ny[0] = nym
        self.dx[0] = dxm
        self.dy[0] = dym
        
        for ig in range(1, nlevel):
            ifg = ig - 1
            self.nx[ig] = int( self.nx[ifg]/2 )
            self.ny[ig] = int( self.ny[ifg]/2 )
            self.dx[ig] = 2 * self.dx[ifg]
            self.dy[ig] = 2 * self.dy[ifg]

    def allocateField( self, vh, fh ):
        # allocate matrix for storage at fine level
        self.eh = np.zeros((self.nx[0]+1, self.ny[0]+1), dtype=np.float64)
        self.rh = np.zeros((self.nx[0]+1, self.ny[0]+1), dtype=np.float64)
        self.vh = vh.copy()
        self.fh = fh.copy()
        
        # allocate matrix for storage at coarse levels
        self.f2h = np.zeros((self.nx[1]+1, self.ny[1]+1), dtype=np.float64)
        self.v2h = np.zeros((self.nx[1]+1, self.ny[1]+1), dtype=np.float64)
        
        self.r2h = np.zeros((self.nx[1]+1, self.ny[1]+1), dtype=np.float64)
        self.e2h = np.zeros((self.nx[1]+1, self.ny[1]+1), dtype=np.float64)
        
        self.f4h = np.zeros((self.nx[2]+1, self.ny[2]+1), dtype=np.float64)
        self.v4h = np.zeros((self.nx[2]+1, self.ny[2]+1), dtype=np.float64)
        
        
    def Multigrid( self ):
        # call relaxation on fine grid and compute the numerical solution
        # for fixed number of iteration
        gauss_seidel_mg(self.nx[0], self.ny[0], self.dx[0], self.dy[0], self.fh, self.vh, self.v1) 
        
        # check for convergence only for finest grid
        # compute the residual and L2 norm
        compute_residual( self.nx[0], self.ny[0], self.dx[0], self.dy[0], self.fh, self.vh, self.rh )
        
        # restrict the residual from fine level to coarse level
        restriction(self.nx[0], self.ny[0], self.nx[1], self.ny[1], self.rh, self.f2h)
        
        # set solution zero on coarse grid
        self.v2h = np.zeros((self.nx[1]+1, self.ny[1]+1), dtype=np.float64)
    
        # solve on the coarsest level and relax V3 times
        gauss_seidel_mg(self.nx[1], self.ny[1], self.dx[1], self.dy[1], self.f2h, self.v2h, self.v3)
        
        #calc r2h
        compute_residual( self.nx[1], self.ny[1], self.dx[1], self.dy[1], self.f2h, self.v2h, self.r2h )
        
        #restric r2h to f4h
        restriction(self.nx[1], self.ny[1], self.nx[2], self.ny[2], self.r2h, self.f4h)
        
        #zero v4h
        self.v4h = np.zeros((self.nx[2]+1, self.ny[2]+1), dtype=np.float64)
        
        #solve v4h
        gauss_seidel_mg(self.nx[2], self.ny[2], self.dx[2], self.dy[2], self.f4h, self.v4h, self.v3)
        
        #interpolate e2h from v4h
        prolongation( self.nx[2], self.ny[2], self.nx[1], self.ny[1], self.v4h, self.e2h )
        
        # correct the solution on fine level
        for j in range(1, self.ny[1]):
            for i in range(1, self.nx[1]):
                self.v2h[i,j] = self.v2h[i,j] + self.e2h[i,j]
        
        # prolongate solution from coarse level to fine level
        prolongation( self.nx[1], self.ny[1], self.nx[0], self.ny[0], self.v2h, self.eh )
        
        # correct the solution on fine level
        for j in range(1, self.ny[0]):
            for i in range(1, self.nx[0]):
                self.vh[i,j] = self.vh[i,j] + self.eh[i,j]
    
        # relax v2 times
        gauss_seidel_mg(self.nx[0], self.ny[0], self.dx[0], self.dy[0], self.fh, self.vh, self.v2)
        
    def init_norm( self ):
        compute_residual(self.nx[0], self.ny[0], self.dx[0], self.dy[0], self.fh, self.vh, self.rh)
        rms = compute_l2norm(self.nx[0], self.ny[0], self.rh)
        return rms

    def compute_norm( self ):
        rms = compute_l2norm(self.nx[0], self.ny[0], self.rh)
        return rms

    def compute_max_residual( self ):
        max_r  = np.max( np.abs(self.rh) )
        return max_r
        
    def Run( self ):
        tolerance = 1.0e-12
        max_iter = 100000
        
        # create output file for L2-norm
        output = open("output.txt", "w");
        output.write("Residual details: \n");    
       
        # create text file for writing residual history
        residual_plot = open("mg_residual.txt", "w")
        rms = self.init_norm()
        init_rms = rms
        
        iter_count = 0
        print(iter_count, " ", init_rms, " ", rms/init_rms)
        
        self.res = []
        self.res.append( rms/init_rms )
            
        # start main iteration loop
        while True:
            self.Multigrid()
            # compute the l2norm of residual
            rms = self.compute_norm()
            
            residual_plot.write("{0} {1} {2}\n".format(iter_count,rms, rms/init_rms));
            print(iter_count, " ", rms, " ", rms/init_rms)
            self.res.append( rms/init_rms )
    
            if (rms/init_rms) <= tolerance or iter_count >= max_iter:
                break
            iter_count += 1
            
        max_r  = self.compute_max_residual()
        output.write("L-2 Norm = {0}\n" .format(rms));
        output.write("Maximum Norm = {0}\n".format(max_r));
        output.write("Iterations =  {0}\n".format(iter_count));
        
        residual_plot.close()
    
def compute_l2norm( nx, ny, r ):
    rms = 0.0
    for j in range(0, ny+1):
        for i in range(0, nx+1):
            rms = rms + r[i,j] * r[i,j]
    rms = np.sqrt( rms / ( ( nx + 1 ) * ( ny + 1 ) ) )
    return rms
    
def compute_residual(nx, ny, dx, dy, f, un, r):
    dx2 = dx * dx
    dy2 = dy * dy
    for j in range(1, ny):
        for i in range(1, nx):
            d2udx2 = (un[i-1,j] - 2*un[i,j] + un[i+1,j])/(dx2)
            d2udy2 = (un[i,j-1] - 2*un[i,j] + un[i,j+1])/(dy2)
            r[i,j] = f[i,j]  - d2udx2 - d2udy2

def gauss_seidel_mg(nx, ny, dx, dy, f, un, niter):
    r = np.zeros( ( nx + 1, ny + 1 ), dtype=np.float64 )
   
    dx2 = dx * dx
    dy2 = dy * dy
    den = -2.0/dx2 - 2.0/dy2
    for n in range( niter ):
        for j in range(1, ny):
            for i in range(1, nx):
                d2udx2 = (un[i-1,j] - 2*un[i,j] + un[i+1,j])/(dx2)
                d2udy2 = (un[i,j-1] - 2*un[i,j] + un[i,j+1])/(dy2)
                r[i,j] = f[i,j] - d2udx2 - d2udy2
                un[i,j] = un[i,j] + r[i,j] / den
                
def restriction(nxf, nyf, nxc, nyc, r, ec):
    for j in range(1, nyc):
        for i in range(1, nxc):
            # grid index for fine grid for the same coarse point
            center = 4.0*r[2*i, 2*j]
            # W, E, N, S with respect to coarse grid point in fine grid
            grid = 2.0*(r[2*i-1, 2*j] + r[2*i+1, 2*j] +
                        r[2*i, 2*j+1] + r[2*i, 2*j-1])
            # NW, NE, SE, SW with respect to coarse grid point in fine grid
            corner = 1.0*(r[2*i-1, 2*j+1] + r[2*i+1, 2*j+1] +
                          r[2*i-1, 2*j-1] + r[2*i+1, 2*j-1])
            # restriction using trapezoidal rule
            ec[i,j] = (center + grid + corner)/16.0

    # restriction for boundary points bottom and top
    #print("nxc,nyc=",nxc,nyc)
    for i in range(0, nxc+1):
        # bottom boundary j = 1
        ec[i, 0] = r[2*i, 0]
        # top boundary j = ny_coarse+1
        ec[i,nyc] = r[2*i, nyf]

    # restriction for boundary poinys left and right
    for j in range(0, nyc+1):
        # left boundary i = 1
        ec[0,j] = r[0,2*j]
        # right boundary nx_coarse+1
        ec[nxc,j] = r[nxf, 2*j]

def prolongation(nxc, nyc, nxf, nyf, unc, ef):
    #print("nxc,nyc,nxf,nyf=",nxc,nyc,nxf,nyf)
    for j in range(0, nyc):
        for i in range(0, nxc):
            # direct injection at center point
            ef[2*i, 2*j] = unc[i,j]
            # east neighnour on fine grid corresponding to coarse grid point
            ef[2*i, 2*j+1] = 0.5*(unc[i,j] + unc[i,j+1])
            # north neighbout on fine grid corresponding to coarse grid point
            ef[2*i+1, 2*j] = 0.5*(unc[i,j] + unc[i+1,j])
            # NE neighbour on fine grid corresponding to coarse grid point
            ef[2*i+1, 2*j+1] = 0.25*(unc[i  ,j] + unc[i  ,j+1] +
                                     unc[i+1,j] + unc[i+1,j+1])

    # update boundary points
    for j in range(0, nyc):
        # left boundary i = 0
        ef[0,2*j+1] = 0.5 * ( unc[0,j] + unc[0,j+1] )
        # right boundary i = nx_fine
        ef[nxf,2*j+1] = 0.5 * ( unc[nxc,j] + unc[nxc,j+1] )
        
    for j in range(0, nyc+1):
        ef[0  ,2*j] = unc[0,j]
        ef[nxf,2*j] = unc[nxc,j]
    
    for i in range(0, nxc):
        #bottom boundary j = 0
        ef[2*i+1,0] = 0.5*(unc[i,0] + unc[i+1,0])
        # top boundary j = ny_fine
        ef[2*i+1,nyf] = 0.5*(unc[i,nyc] + unc[i+1,nyc])

    for i in range(0, nxc+1):
        ef[2*i,0] = unc[i,0]
        ef[2*i,nyf] = unc[i,nyc]

#n=128
#n=4
n=64
print("n=",n)
nx = n
ny = n

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
    
ue = np.zeros( ( nx + 1, ny + 1 ) )
un = np.zeros( ( nx + 1, ny + 1 ) )
f = np.zeros( ( nx + 1, ny + 1 ) )
r = np.zeros( ( nx + 1, ny + 1 ) )

c1 = np.power(1.0/16.0,2)
c2 = - 2.0 * np.pi * np.pi

for j in range(0, ny):
    for i in range(0, nx):
      ue[i,j] = np.sin(2.0*np.pi*x[i])*np.sin(2.0*np.pi*y[j]) + \
                c1*np.sin(16.0*np.pi*x[i])*np.sin(16.0*np.pi*y[j])
      
      f[i,j] = 4.0*c2*np.sin(2.0*np.pi*x[i])*np.sin(2.0*np.pi*y[j]) + \
               c2*np.sin(16.0*np.pi*x[i])*np.sin(16.0*np.pi*y[j])
      
      un[i,j] = 0.0
    
nlevel = 3
mg = Mgdata()
mg.allocate( nx, ny, dx, dy, nlevel )
mg.allocateField( un, f )
mg.Run()

fig = plt.figure("OneFLOW-CFD Solver+Restriction", figsize=(10, 6), dpi=100)

ax = fig.add_subplot(1, 2, 1)
ax.axis('off')

cntr = ax.contourf(x, y, mg.vh, levels=20, cmap="jet")
plt.colorbar(cntr,ax=ax)
ax.title.set_text("Exact solution")

ax = fig.add_subplot(1, 2, 2)
nCycle = len( mg.res )
ii = np.linspace(0, nCycle-1, nCycle)
ax.set_yscale("log", base=10)
ax.plot(ii,mg.res,color='black',linewidth=2,linestyle='solid',marker='none',markerfacecolor='none',label='res')

plt.tight_layout()
plt.show()
