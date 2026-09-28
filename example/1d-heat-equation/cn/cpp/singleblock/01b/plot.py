import numpy as np
import matplotlib.pyplot as plt
import csv

x_list = []
ue_list = []
un_list = []
uerror_list = []

with open('field_final.csv', newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for row in readCSV:
        #print(row)
        #print(row[0])
        #print(row[0], row[1], row[2], row[3])
        if ( icount != 0 ) :
            x_list.append(row[0])
            ue_list.append(row[1])
            un_list.append(row[2])
            uerror_list.append(row[3])
        icount += 1

ni = icount - 1
print("ni=",ni)

x = np.zeros( ni )
ue = np.zeros( ni )
un = np.zeros( ni )
uerror = np.zeros( ni )

for i in range(0, ni):
    x[i] = float(x_list[i])
    ue[i] = float(ue_list[i])
    un[i] = float(un_list[i])
    uerror[i] = float(uerror_list[i])
    
plt.figure("OneFLOW-CFD Solver", figsize=(10, 4), dpi=100)
plt.subplot(1, 2, 1)
plt.plot(x, ue, "k-", linewidth=1.0, label="Exact solution")
plt.scatter(x, un, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="CN solution")
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
    

