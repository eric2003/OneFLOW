import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
def compute_l2norm(nx,r):
    rms = 0.0
    for i in range(1, nx):
        rms += r[i] * r[i]
    rms = np.sqrt( rms / ( ( nx - 1 ) ) )
    return rms
    
def compute_max_error( nx, u_error ):
    val_max = -1;
    ipos = -1;
    for i in range(1, nx):
        if ( val_max < np.abs( u_error[ i ] ) ):
            ipos = i;
            val_max = np.abs( u_error[ i ] )
    print( " ipos = ", ipos )
    return val_max;


nt = 400
t  = 1.0

filename = 'field_final'+str(nt)+'.csv'

nvar = len(sys.argv)
print('nvar=',nvar)
print('sys.argv=',sys.argv)

labelname = "FTCS solution"
if nvar >= 2:
    scheme = sys.argv[1]
    print('scheme=',scheme)
    labelname = scheme + ' solution'
    
print("labelname=",labelname)    

with open(filename, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for row in readCSV:
        icount += 1

ni = icount
print("ni=",ni)

x = np.zeros( (ni) )
u = np.zeros( (ni) )
ue= np.zeros( (ni) )

with open(filename, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    i = 0
    for row in readCSV:
        x[i] = float(row[0])
        u[i] = float(row[1])
        i += 1
#sort
sorted_indices = np.argsort(x)
xt=x[sorted_indices]
ut=u[sorted_indices]
x=xt
u=ut

for i in range(ni):
    ue[i] = - np.exp(-t) * np.sin( np.pi * x[i] ) # theory solution 
    
uerror = u - ue

nx = ni - 1

my_max_error = compute_max_error( nx, uerror )
print("my_max_error = {0:.15f}".format(my_max_error))

rms_error = compute_l2norm(nx,uerror)
max_error = np.max( np.abs(uerror) )

# create output file for L2-norm
output = open("output.txt", "w");
output.write("Error details: \n");
output.write("L-2 Norm = {0}\n" .format(str(rms_error)));
output.write("Maximum Norm = {0}\n".format(str(max_error)));
output.close()
  
plt.figure("OneFLOW-CFD Solver", figsize=(10, 4), dpi=100)
plt.subplot(1, 2, 1)
plt.plot(x, ue, "k-", linewidth=1.0, label="Exact solution")
#plt.scatter(x, u, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="ICP solution")
plt.scatter(x, u, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=labelname)
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Solution field")
plt.legend()
plt.tight_layout()

plt.subplot(1, 2, 2)
plt.scatter(x, np.abs(uerror), facecolor="none", edgecolor="green", s=20, linewidths=0.5)
plt.ylabel(r"$\epsilon$")
plt.xlabel("$x$")
plt.title("Discretization error")
plt.tight_layout()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))

plt.show();    
    

