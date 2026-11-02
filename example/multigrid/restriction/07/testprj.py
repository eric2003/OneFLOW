import matplotlib.pyplot as plt
import numpy as np

class Mgdata:
    def __init__(self):
        self.nlevel = 1
    
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
            
def MultigridCycle(nx, ny, dx, dy, un, f, r, output):
    tolerance = 1.0e-12
    max_iter = 100000
    nlevel = 2
    mg = Mgdata()
    mg.allocate( nx, ny, dx, dy, nlevel )
    #lnx, lny, ldx, ldy = allocate( nlevel, nx, ny, dx, dy )
    
    # allocate matrix for storage at fine level
    # residual at fine level is already defined at global level
    prol_fine = np.zeros((mg.nx[0]+1, mg.ny[0]+1), dtype=np.float64)
    
    # allocate matrix for storage at coarse levels
    fc  = np.zeros((mg.nx[1]+1, mg.ny[1]+1), dtype=np.float64)
    unc = np.zeros((mg.nx[1]+1, mg.ny[1]+1), dtype=np.float64)
   
    v1 = 2 # relaxation
    v2 = 2 # prolongation
    v3 = 2 # coarsest level    
    
    # create text file for writing residual history
    residual_plot = open("mg_residual.txt", "w")

    compute_residual(nx, ny, dx, dy, f, un, r)

    rms = compute_l2norm(nx, ny, r)

    init_rms = rms
    iter_count = 0
    print(iter_count, " ", init_rms, " ", rms/init_rms)
    
    res = []
    res.append( rms/init_rms )
        
    # start main iteration loop
    while True:
        # call relaxation on fine grid and compute the numerical solution
        # for fixed number of iteration
        gauss_seidel_mg(mg.nx[0], mg.ny[0], mg.dx[0], mg.dy[0], f, un, v1) 
        
        # check for convergence only for finest grid
        # compute the residual and L2 norm
        compute_residual( mg.nx[0], mg.ny[0], mg.dx[0], mg.dy[0], f, un, r )
        
        # compute the l2norm of residual
        rms = compute_l2norm(mg.nx[0], mg.ny[0], r)
        
        residual_plot.write("{0} {1} {2}\n".format(iter_count,rms, rms/init_rms));
        print(iter_count, " ", rms, " ", rms/init_rms)
        res.append( rms/init_rms )

        if (rms/init_rms) <= tolerance or iter_count >= max_iter:
            break
        iter_count += 1
        
        # restrict the residual from fine level to coarse level
        restriction(mg.nx[0], mg.ny[0], mg.nx[1], mg.ny[1], r, fc)
        
        # set solution zero on coarse grid
        unc = np.zeros((mg.nx[1]+1, mg.ny[1]+1), dtype=np.float64)

        # solve on the coarsest level and relax V3 times
        gauss_seidel_mg(mg.nx[1], mg.ny[1], mg.dx[1], mg.dy[1], fc, unc, v3)
        
        # prolongate solution from coarse level to fine level
        prolongation( mg.nx[1], mg.ny[1], mg.nx[0], mg.ny[0], unc, prol_fine )
        
        # correct the solution on fine level
        for j in range(1, mg.ny[0]):
            for i in range(1, mg.nx[0]):
                un[i,j] = un[i,j] + prol_fine[i,j]

        # relax v2 times
        gauss_seidel_mg(mg.nx[0], mg.ny[0], mg.dx[0], mg.dy[0], f, un, v2) 
        
    max_r  = np.max( np.abs(r) )  
    output.write("L-2 Norm = {0}\n" .format(rms));
    output.write("Maximum Norm = {0}\n".format(max_r));
    output.write("Iterations =  {0}\n".format(iter_count));
    
    residual_plot.close()
    return res

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
    
# create output file for L2-norm
output = open("output.txt", "w");
output.write("Residual details: \n");

res = MultigridCycle(nx, ny, dx, dy, un, f, r, output )

fig = plt.figure("OneFLOW-CFD Solver+Restriction", figsize=(10, 6), dpi=100)

ax = fig.add_subplot(1, 2, 1)
ax.axis('off')

cntr = ax.contourf(x, y, un, levels=20, cmap="jet")
plt.colorbar(cntr,ax=ax)
ax.title.set_text("Exact solution")

ax = fig.add_subplot(1, 2, 2)
nCycle = len( res )
ii = np.linspace(0, nCycle-1, nCycle)
ax.set_yscale("log", base=10)
ax.plot(ii,res,color='black',linewidth=2,linestyle='solid',marker='none',markerfacecolor='none',label='res')

plt.tight_layout()
plt.show()
