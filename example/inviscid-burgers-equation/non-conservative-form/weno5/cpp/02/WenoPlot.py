import numpy as np
import matplotlib.pyplot as plt
import csv

with open('field_final.csv', newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    headers = next(readCSV)
    num_columns = len(headers) 
    print(f'num_columns = {num_columns}')
    icount = 0
    for row in readCSV:
        icount += 1

ni = icount + 1
print("ni=",ni)

ns = num_columns - 2

u = np.zeros( (ni, ns + 1 ) )
x = np.zeros( (ni) )

with open('field_final.csv', newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    i = 0
    for row in readCSV:
        #print(f"row={row}\n")
        x[i] = float(row[0])
        for j in range(1,num_columns):
            u[i][j-1] = float(row[j])
        i += 1


print("u.shape=",u.shape)
n1 = u.shape[0]
n2 = u.shape[1]
print(f"n1={n1},n2={n2}")

#x = np.linspace(0,1, num=ni) 

tm = 0.25

plt.figure("OneFLOW-CFD Solver", figsize=(6, 4), dpi=100)
for k in range(0, ns+1):
    plt.plot(x, u[:,k], linewidth=1.0, label="t="+format(tm*k/ns, ".4f"))
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Inviscid Burgers Equation: Non-Conservative Form-WENO-5 Scheme")
plt.legend(loc='upper right', fontsize='6')
plt.show()    