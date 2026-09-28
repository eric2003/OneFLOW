import numpy as np
import matplotlib.pyplot as plt
import csv

with open('field_final0.csv', newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for row in readCSV:
        icount += 1

ni = icount
print("ni=",ni)

ns = 10

u = np.zeros( (ni, ns + 1 ) )
x = np.zeros( (ni) )

for j in range(ns+1):
    filename = 'field_final'+str((j)*250)+'.csv'
    print('filename=',filename)
    with open(filename, newline='') as csvfile:
        readCSV = csv.reader(csvfile, delimiter= ' ')
        i = 0
        for row in readCSV:
            x[i] = float(row[0])
            u[i][j] = float(row[1])
            i += 1
    
print("u.shape=",u.shape)
n1 = u.shape[0]
n2 = u.shape[1]
print(f"n1={n1},n2={n2}")
#exit()
#x = np.linspace(0,1, num=ni) 

#sort
sorted_indices = np.argsort(x)
xt=x[sorted_indices]

for k in range(ns+1):
    ut = u[sorted_indices,k]
    u[:,k] = ut[:]
tm = 0.25

plt.figure("OneFLOW-CFD Solver", figsize=(6, 4), dpi=100)
for k in range(0, ns+1):
    plt.plot(xt, u[:,k], linewidth=1.0, label="t="+format(tm*k/ns, ".4f"))
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Inviscid Burgers Equation: Non-Conservative Form-WENO-5 Scheme")
plt.legend(loc='upper right', fontsize='6')
plt.show()    