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
# Calculate right hand term
#-----------------------------------------------------------------------------#
def rhs( nx, dx, u, r, alpha ):
    for i in range(1, nx):
        r[i] = alpha * ( u[i+1] - 2.0*u[i] + u[i-1] )/( dx * dx )
    return
    
#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#-----------------------------------------------------------------------------#    
def numerical( nx, nt, dx, dt, x, u, alpha ):
    un = np.zeros( (nx+1) ) # numerical solution at every time step
    ut = np.zeros( (nx+1) ) # temporary array during RK3 integration
    
    r = np.zeros( (nx+1) )

    for i in range(0, nx+1):
        un[i] = - np.sin( np.pi * x[i] )
        u[i,0] = un[i] # store solution at t=0
    
    # dirichlet boundary condition
    un[0 ] = 0.0
    un[nx] = 0.0

    # dirichlet boundary condition for temporary array
    ut[0 ] = 0.0
    ut[nx] = 0.0

    for k in range(1, nt+1):
        rhs( nx, dx, un, r, alpha )

        for i in range(1, nx):
            ut[i] = un[i] + dt * r[i]

        rhs( nx, dx, ut, r, alpha )

        for i in range(1, nx):
            ut[i] = 0.75 * un[i] + 0.25 * ut[i] + 0.25 * dt * r[i]

        rhs( nx, dx, ut, r, alpha )

        for i in range(1, nx):
            un[i] = ( 1.0 / 3.0 ) * un[i] + ( 2.0 / 3.0 ) * ut[i] + ( 2.0 / 3.0 ) * dt * r[i]

        u[:,k] = un[:]
    return

x_l = -1.0
x_r = 1.0
dx = 0.025
nx = int( ( x_r - x_l ) / dx )

print( " nx = ", nx )

dt = 0.0025
t  = 1.0
nt = int( t / dt )
print( " nt = ", nt )

alpha = 1 / ( np.pi**2 )
print( " alpha = ", alpha )

x = np.zeros( (nx+1) )  
u_e = np.zeros( (nx+1) )  
u = np.zeros( (nx+1, nt+1 ) )  

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)  # location of each grid point
    u_e[i] = - np.exp(-t) * np.sin( np.pi * x[i] ) # theory solution
    
numerical( nx, nt, dx, dt, x, u, alpha )    

# compute L2 norm of the error
u_error = u[:,nt] - u_e
rms_error = compute_l2norm(nx,u_error)
max_error = np.max( np.abs(u_error) )

# create output file for L2-norm
output = open("output.txt", "w");
output.write("Error details: \n");
output.write("L-2 Norm = {0}\n" .format(str(rms_error)));
output.write("Maximum Norm = {0}\n".format(str(max_error)));
output.close()

plt.figure("OneFLOW-CFD Solver", figsize=(10, 4), dpi=100)
plt.subplot(1, 2, 1)
plt.plot(x, u_e, "k-", linewidth=1.0, label="Exact solution")
plt.scatter(x, u[:,nt], facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="Runge-Kutta 3 scheme")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Solution field")
plt.legend()
plt.tight_layout()

plt.subplot(1, 2, 2)
plt.plot(x, np.abs(u_error), color="green", marker="o", linewidth=1.0,label="Error")
plt.ylabel(r"$\epsilon$")
plt.xlabel("$x$")
plt.title("Discretization error")
plt.legend()
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))

plt.show();