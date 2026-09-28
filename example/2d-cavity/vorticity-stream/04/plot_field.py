import numpy as np
import matplotlib.pyplot as plt

x_list = []
y_list = []
wn_list = []
sn_list = []

with open('field_final.txt', 'r') as f:
    for index, line in enumerate(f):
        words = line.strip().split()
        if index == 0:
            nx = int(words[0])
            ny = int(words[1])
        else:
            x_list.append( float(words[0]) )
            y_list.append( float(words[1]) )
            wn_list.append( float(words[2]) )
            sn_list.append( float(words[3]) )
            
x  = np.zeros( ( nx + 1, ny + 1 ) )
y  = np.zeros( ( nx + 1, ny + 1 ) )
wn = np.zeros( ( nx + 1, ny + 1 ) )
sn = np.zeros( ( nx + 1, ny + 1 ) )

id = 0        
for j in range(0, ny+1):
    for i in range(0, nx+1):
        x[i,j] = x_list[id]
        y[i,j] = y_list[id]
        wn[i,j] = wn_list[id]
        sn[i,j] = sn_list[id]
        id += 1
        
# Creating 2-D grid of features
#[X, Y] = np.meshgrid(x, y, indexing='ij')

fig = plt.figure("OneFLOW-CFD ", figsize=(14, 6), dpi=100)
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

cs = ax1.contourf(x, y, wn, levels=60, cmap="jet")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_title("Vorticity field")
fig.colorbar(cs, ax = ax1)

cs = ax2.contourf(x, y, sn, levels=60, cmap="jet")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_title("Streamfunction")
fig.colorbar(cs, ax = ax2)
fig.tight_layout()
fig.savefig("ldc_contour.pdf")
plt.show()

