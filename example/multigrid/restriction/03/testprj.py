import matplotlib.pyplot as plt
import numpy as np

def plot_point( xp, yp, ax ):
    x = [xp]
    y = [yp]
    ax.plot(x, y, marker="o", markersize=8, markeredgecolor="blue",
    markerfacecolor="blue")
    
def FieldNorm(vh): 
    return np.max(np.abs(vh))
    
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

    
def gauss_seidel(un, f, r, nx, ny, dx, dy, niter):
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
            
def MultigridCycle(nx, ny, dx, dy, un, f, r, output):
    tolerance = 1.0e-10
    max_iter = 100000
    #allocate memory for grid size at different levels
    lnx = np.zeros(2, dtype=np.int32 )
    lny = np.zeros(2, dtype=np.int32 )
    ldx = np.zeros(2, dtype=np.float64)
    ldy = np.zeros(2, dtype=np.float64)
    
    lnx[0] = nx
    lny[0] = ny
    lnx[1] = int( lnx[0]/2 )
    lny[1] = int( lny[0]/2 )
    ldx[0] = dx
    ldy[0] = dy
    ldx[1] = 2 * ldx[0]
    ldy[1] = 2 * ldy[0]
    print("lnx=\n",lnx)
    print("lny=\n",lny)
    print("ldx=\n",ldx)
    print("ldy=\n",ldy)
    #exit()
    
    # allocate matrix for storage at fine level
    # residual at fine level is already defined at global level
    prol_fine = np.zeros((lnx[0]+1, lny[0]+1), dtype=np.float64)
    
    # allocate matrix for storage at coarse levels
    fc  = np.zeros((lnx[1]+1, lny[1]+1), dtype=np.float64)
    unc = np.zeros((lnx[1]+1, lny[1]+1), dtype=np.float64)
    
    niter = 3
    
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
        gauss_seidel_mg(lnx[0], lny[0], dx, dy, f, un, niter) 
        
        # check for convergence only for finest grid
        # compute the residual and L2 norm
        compute_residual( nx, ny, dx, dy, f, un, r )
        
        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, r)
        
        residual_plot.write("{0} {1} {2}\n".format(iter_count,rms, rms/init_rms));
        print(iter_count, " ", rms, " ", rms/init_rms)
        res.append( rms/init_rms )

        if (rms/init_rms) <= tolerance or iter_count >= max_iter:
            break
        iter_count += 1
        
        # relax v2 times
        gauss_seidel_mg(lnx[0], lny[0], dx, dy, f, un, niter) 
        
        
    max_r  = np.max( np.abs(r) )  
    output.write("L-2 Norm = {0}\n" .format(rms));
    output.write("Maximum Norm = {0}\n".format(max_r));
    output.write("Iterations =  {0}\n".format(iter_count));
    
    residual_plot.close()
    return res

#n=4
n=16
nx = n
ny = n
h=1/n
print("n=",n)
print("h=",h)
dx = h
dy = h

x = np.zeros(n+1)
y = np.zeros(n+1)
xmin = 0
xmax = 1
ymin = 0
ymax = 1
x[0] = xmin
for i in range(n):
    #print("i=",i,"n=",n)
    x[i+1] = xmin + ( i + 1 ) * dx

for j in range(n):
    y[j+1] = ymin + ( j + 1 ) * dy
    
ue = np.zeros( ( n + 1, n + 1 ) )
un = np.zeros( ( n + 1, n + 1 ) )
f = np.zeros( ( n + 1, n + 1 ) )
r = np.zeros( ( n + 1, n + 1 ) )

for j in range(1, n):
    for i in range(0, n+1):
        y16m = np.sin( i * 16 * np.pi / n )
        y40m = np.sin( i * 40 * np.pi / n )
        #y16m = np.sin( i * np.pi / n )
        #y40m = np.sin( i * np.pi / n )
        ue[i,j] = 0.5 * ( y16m + y40m )
        un[i,j] = ue[i,j]
        
# create output file for L2-norm
output = open("output.txt", "w");
output.write("Residual details: \n");

res = MultigridCycle(nx, ny, dx, dy, un, f, r, output )

#print("res=\n",res)
    
fig = plt.figure("OneFLOW-CFD Solver+Restriction", figsize=(6, 8), dpi=100)

ax = fig.add_subplot(2, 2, 1)

ax.title.set_text("Line Value")

for j in range(n+1):
    #print("j=",j,"un=\n",un[:,j])
    ax.plot(x,un[:,j],color='black',linewidth=1)

ax = fig.add_subplot(2, 2, 2)
ax.axis('off')

for j in range(n+1):
    ax.hlines(y = y[j], xmin = x[0], xmax = x[-1], colors="k" )

for i in range(n+1):
    ax.vlines(x = x[i], ymin = y[0], ymax = y[-1], colors="k")
    

cntr = ax.contourf(x, y, un, levels=20, cmap="jet")
plt.colorbar(cntr,ax=ax)
ax.title.set_text("Exact solution")

ax = fig.add_subplot(2, 2, 3)
nCycle = len( res )
ii = np.linspace(0, nCycle-1, nCycle)
ax.set_yscale("log", base=10)
ax.plot(ii,res,color='black',linewidth=2,linestyle='solid',marker='none',markerfacecolor='none',label='res')

plt.tight_layout()
plt.show()
