import numpy as np
import matplotlib.pyplot as plt
import csv

x_list = []
u_list = []

with open('solution_total.plt', newline='') as csvfile:
    #readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for line in csvfile:
        # 去除首尾空格，按连续空格分割
        row = line.strip().split()    
        #print("row=",row)
        x_list.append(row[0])
        u_list.append(row[1])
        icount += 1

ni = icount
print("ni=",ni)

#print("x_list=",x_list)

x = np.zeros( ni )
u = np.zeros( ni )

for i in range(0, ni):
    x[i] = float(x_list[i])
    u[i] = float(u_list[i])
    
plt.figure("OneFLOW-CFD Solver", figsize=(8, 6), dpi=100)
plt.plot(x, u, "b-", linewidth=1.0, marker = 's', markerfacecolor='None', label="Eno2")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("burgers T=1")
plt.legend()
plt.tight_layout()


plt.show();    
    

