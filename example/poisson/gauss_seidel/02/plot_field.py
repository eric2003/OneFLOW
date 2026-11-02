import numpy as np
import matplotlib.pyplot as plt

x_list = []
y_list = []
f_list = []
un_list = []
ue_list = []

with open('field_final.txt', 'r') as f:
    for index, line in enumerate(f):
        words = line.strip().split()
        if index == 0:
            nx = int(words[0])
            ny = int(words[1])
        else:
            x_list.append( float(words[0]) )
            y_list.append( float(words[1]) )
            f_list.append( float(words[2]) )
            un_list.append( float(words[3]) )
            ue_list.append( float(words[4]) )
            
x = np.zeros( ( nx + 1, ny + 1 ) )
y = np.zeros( ( nx + 1, ny + 1 ) )
f  = np.zeros( ( nx + 1, ny + 1 ) )
un = np.zeros( ( nx + 1, ny + 1 ) )
ue = np.zeros( ( nx + 1, ny + 1 ) )

id = 0        
for j in range(0, ny+1):
    for i in range(0, nx+1):
        x[i,j] = x_list[id]
        y[i,j] = y_list[id]
        f[i,j] = f_list[id]
        un[i,j] = un_list[id]
        ue[i,j] = ue_list[id]
        id += 1


fig = plt.figure("OneFLOW-CFD Solver+2D Poisson Equation+Gauss-Seidel", figsize=(10, 4), dpi=100)

ax1 = fig.add_subplot(1,2,1)
cs1 = ax1.contourf( x, y, ue, levels=20, cmap="jet",vmin=-1,vmax=1)
fig.colorbar(cs1, ax = ax1)
ax1.set_title("Exact solution")

ax2 = fig.add_subplot(1,2,2)
cs2 = ax2.contourf( x, y, un, levels=20, cmap="jet",vmin=-1,vmax=1)
fig.colorbar(cs2, ax = ax2)
ax2.set_title("Numerical solution")

fig.tight_layout()
fig.savefig("gs_contour.pdf")
plt.show()
