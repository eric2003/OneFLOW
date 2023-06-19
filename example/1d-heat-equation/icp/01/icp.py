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

a = np.zeros(nx+1)
b = np.zeros(nx+1)
c = np.zeros(nx+1)
d = np.zeros(nx+1)

half = 1.0/2.0

r = half * alpha * dt / ( dx * dx )

for k in range(1, nt+1):
    for i in range(0, nx + 1):
        a[i] = 1.0 /12.0 - r
        b[i] = 10.0/12.0 + 2.0 * r
        c[i] = 1.0 /12.0 - r
        
    a[0] = 0
    b[0] = 1
    c[0] = 0
    
    a[nx] = 0
    b[nx] = 1
    c[nx] = 0
    
    for i in range(1, nx):
        aa = 1.0 /12.0 + r
        bb = 10.0/12.0 - 2.0 * r
        cc = 1.0 /12.0 + r
        d[i] = aa * u[i-1,k-1] + bb * u[i,k-1] + cc * u[i+1,k-1]
        
    d[0] = 0
    d[nx] = 0
        
    value = thomas_algorithm( a, b, c, d, nx + 1 )
    
    u[:,k] = value
        

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
plt.scatter(x, u[:,nt], facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="Implicit Compact Pade Scheme")
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
#plt.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))
plt.ticklabel_format(axis='y', style='sci', scilimits=(-8,-8))

plt.show();