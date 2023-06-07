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
un = np.zeros( (nt+1, nx+1) )  
error = np.zeros( (nt+1, nx+1) )  

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)  # location of each grid point
    un[0,i] = - np.sin( np.pi * x[i] ) # initial condition @ t=0
    u_e[i] = - np.exp(-t) * np.sin( np.pi * x[i] ) # theory solution

un[0,0 ] = 0.0
un[0,nx] = 0.0

beta = alpha * dt / ( dx**2 )

for k in range(1, nt+1):
    for i in range(1, nx):
        un[k,i] = un[k-1,i] + beta * ( un[k-1,i+1] - 2.0 * un[k-1,i] + un[k-1,i-1] )
    un[k,0 ] = 0.0 # boundary condition at x = -1
    un[k,nx] = 0.0 # boundary condition at x = 1

# compute L2 norm of the error
u_error = un[nt,:] - u_e
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
plt.scatter(x, un[nt, :], facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="FTCS solution")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Solution field")
plt.legend()
plt.tight_layout()

plt.subplot(1, 2, 2)
plt.scatter(x, np.abs(u_error), facecolor="none", edgecolor="green", s=20, linewidths=0.5)
plt.ylabel(r"$\epsilon$")
plt.xlabel("$x$")
plt.title("Discretization error")
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))

plt.show();