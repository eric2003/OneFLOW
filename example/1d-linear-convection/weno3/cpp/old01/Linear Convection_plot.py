import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

def set_fuction(nx, u, x, ct) :
    for i in range(nx):
        xm = x[i]
        if 0.5 <= xm - ct <= 1:
            u[i] = 2
        else:
            u[i] = 1

filename = 'field_final.csv'

nvar = len(sys.argv)
print('nvar=',nvar)
print('sys.argv=',sys.argv)

labelname = "FTBS solution"
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
utheory= np.zeros( (ni) )

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

c = 1
total_t = 0.625
# theory solution 
set_fuction(ni, utheory, x, c * total_t )
  
plt.figure("OneFLOW-CFD Solver")
plt.plot(x, utheory, "k-", linewidth=1.0, label="Exact solution")
plt.scatter(x, u, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label=labelname)
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Solution field")
plt.legend()
plt.tight_layout()
plt.show();

