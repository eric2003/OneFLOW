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
u = np.zeros( (nx+1, nt+1 ) )  

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)  # location of each grid point
    u[i,0] = - np.sin( np.pi * x[i] ) # store solution at t=0   
    u_e[i] = - np.exp(-t) * np.sin( np.pi * x[i] ) # theory solution
    
beta = alpha * dt / ( dx**2 )
print( " beta = ", beta )

#Third-order Runge-Kutta version3(Ralston formula):
    
un = np.zeros(nx+1)
u1 = np.zeros(nx+1)
u2 = np.zeros(nx+1)
k1 = np.zeros(nx+1)
k2 = np.zeros(nx+1)
k3 = np.zeros(nx+1)

half = 1.0/2.0

for k in range(1, nt+1):
    for i in range(1, nx):
        un[i] = u[i,k-1]
        
    for i in range(1, nx):
        k1[i] = beta * ( un[i+1] - 2.0*un[i] + un[i-1] )    
        
    for i in range(1, nx):
        u1[i] = un[i] + ( 1.0 / 2.0) * k1[i]
        
    for i in range(1, nx):
        k2[i] = beta * ( u1[i+1] - 2.0*u1[i] + u1[i-1] )           

    for i in range(1, nx):
        u2[i] = un[i] + ( 3.0 / 4.0 ) * k2[i]

    for i in range(1, nx):
        k3[i] = beta * ( u2[i+1] - 2.0*u2[i] + u2[i-1] )           
        
    for i in range(1, nx):
        u[i,k] = un[i] + ( 1.0 / 9.0 ) * ( 2.0 * k1[i] + 3.0 * k2[i] + 4.0 * k3[i] )

    u[0 ,k] = 0
    u[nx,k] = 0

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