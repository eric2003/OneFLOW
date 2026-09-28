import numpy as np
import matplotlib.pyplot as plt
import time

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
            
def calc_ap( u, ap, dx, dy, nx, ny ):
    dx2 = dx * dx
    dy2 = dy * dy
    
    for j in range(1, ny):
        for i in range(1, nx):
            d2udx2 = ( u[i+1,j] - 2*u[i,j] + u[i-1,j] ) / dx2
            d2udy2 = ( u[i,j+1] - 2*u[i,j] + u[i,j-1] ) / dy2
            ap[i,j] = d2udx2 + d2udy2
    return
    
def my_dot_product( a, b, nx, ny ):
    sum = 0.0
    for j in range(1, ny):
        for i in range(1, nx):
            sum += a[i,j] * b[i,j]
    return sum    

def conjugate_gradient(dx, dy, nx, ny, r, f, u, rms, init_rms, max_iter, tolerance, tiny, output):
    # create text file for writing residual history
    residual_plot = open("cg_residual.txt", "w")
    count = 0.0

    compute_residual(nx, ny, dx, dy, f, u, r)
    rms = compute_l2norm(nx, ny, r)
    
    init_rms = rms
    iter_count = 0
    print(iter_count, " ", init_rms, " ", rms/init_rms)
    
    # allocate the matric for direction and set the initial direction (conjugate vector)
    p = np.zeros( ( nx+1, ny+1) )
    del_p = np.zeros( ( nx+1, ny+1) )
    
    # asssign conjugate vector to initial residual
    for j in range(0, ny):
        for i in range(0, nx):
            p[i,j] = r[i,j]

    while iter_count < max_iter:
        # calculate A*p
        calc_ap( p, del_p, dx, dy, nx, ny )
        
        aa = my_dot_product( r, r, nx, ny )
        bb = my_dot_product( p, del_p, nx, ny )
        
        alpha = aa / ( bb + tiny )

        # update the numerical solution by adding some component of conjugate vector
        for j in range(1, ny):
            for i in range(1, nx):
                u[i,j] = u[i,j] + alpha * p[i,j]
            
        # update the residual by removing some component of previous residual
        for j in range(1, ny):
            for i in range(1, nx):
                r[i,j] = r[i,j] - alpha * del_p[i,j]
                
        # bb = <r(k),r(k)> = aa (calculated in previous loop)
        bb = aa
        aa = my_dot_product( r, r, nx, ny )
        
        # beta = <r(k+1), r(k+1)>/<r(k),r(k)>
        beta = aa / ( bb + tiny )
       
        # update the conjugate vector
        for j in range(0, ny):
            for i in range(0, nx):
                p[i,j] = r[i,j] + beta * p[i,j]
                
        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, r)
        
        residual_plot.write("{0} {1} {2}\n".format(iter_count,rms, rms/init_rms));

        count = iter_count

        print(iter_count, " ", rms, " ", rms/init_rms)
        
        iter_count += 1

        if ( rms/init_rms ) <= tolerance:
            break
        
   
    

nx = 512
ny = 512

tolerance = 1.0e-10
max_iter = 20*100000
tiny = 1.0e-16

print("max_iter",max_iter)

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

dx = (x_r-x_l)/nx
dy = (y_t-y_b)/ny

x = np.zeros( nx + 1 )
y = np.zeros( ny + 1 )
ue = np.zeros( ( nx + 1, ny + 1 ) )
f  = np.zeros( ( nx + 1, ny + 1 ) )
un = np.zeros( ( nx + 1, ny + 1 ) )

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)

for j in range(0, ny+1):
    y[j] = y_b + dy*(j)

# given exact solution
for j in range(0, ny+1):
    for i in range(0, nx+1):
        ue[i,j] = ( x[i] * x[i] - 1.0 ) * ( y[j] * y[j] - 1.0 )
        f[i,j] = -2.0 * ( 2.0 - x[i] * x[i] - y[j] * y[j] )
        un[i,j] = 0.0
        
un[:,0 ] = ue[:,0]
un[:,ny] = ue[:,ny]

un[0 ,:] = ue[0,:]
un[nx,:] = ue[nx,:]

r = np.zeros( ( nx + 1, ny + 1 ) )
init_rms = 0.0
rms = 0.0

# create output file for L2-norm
output = open("output.txt", "w");
output.write("Residual details: \n");
       
time_1 = time.time()

#gauss_seidel(dx, dy, nx, ny, r, f, un, rms, init_rms, max_iter, tolerance, output)
conjugate_gradient(dx, dy, nx, ny, r, f, un, rms, init_rms, max_iter, tolerance, tiny, output)
                      
                      
time_2 = time.time()
t = time_2 - time_1

uerror = np.zeros( ( nx + 1, ny + 1 ) )
rms_error = 0.0

uerror = un - ue

rms_error = compute_l2norm(nx, ny, uerror)
max_error = np.max( np.abs(uerror) )

print("Error details:");
print("L-2 Norm = ", rms_error);
print("Maximum Norm = ", max_error);
print("CPU Time = ", t);

#output = open("output.txt", "w");
output.write("Error details: \n");
output.write("L-2 Norm = {0}\n" .format(rms_error));
output.write("Maximum Norm = {0}\n".format(max_error));
output.write("CPU Time = {0}\n".format(t));
output.close()

titles = ("Exact solution", "Numerical solution")

plt.figure("OneFLOW-CFD Solver+2D Poisson Equation+Conjugate Gradient", figsize=(10, 4), dpi=100)
for i, result in enumerate([ue, un]):
    plt.subplot(1, 2, i+1)
    plt.contourf( x, y, result, levels=20, cmap="jet")
    plt.colorbar()
    plt.title(titles[i])
plt.tight_layout()
plt.show()
