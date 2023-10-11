import numpy as np
import matplotlib.pyplot as plt


def read_field(filename):
    x_list = []
    y_list = []
    wn_list = []
    
    nx = 128
    ny = 128
    
    with open(filename, 'r') as f:
        for index, line in enumerate(f):
            words = line.strip().split()
            x_list.append( float(words[0]) )
            y_list.append( float(words[1]) )
            wn_list.append( float(words[2]) )

                
    x  = np.zeros( ( nx + 1, ny + 1 ) )
    y  = np.zeros( ( nx + 1, ny + 1 ) )
    wn = np.zeros( ( nx + 1, ny + 1 ) )
    
    id = 0        
    for j in range(0, ny+1):
        for i in range(0, nx+1):
            x[i,j] = x_list[id]
            y[i,j] = y_list[id]
            wn[i,j] = wn_list[id]
            id += 1
    return nx,ny,x,y,wn
    
nx,ny,x,y,wn0  = read_field("vm0.txt")
nx,ny,x,y,wn10 = read_field("vm5.txt")
nx,ny,x,y,wn20 = read_field("vm10.txt")
        
# Creating 2-D grid of features
#[X, Y] = np.meshgrid(x, y, indexing='ij')

fig = plt.figure("OneFLOW-CFD ", figsize=(16, 6), dpi=100)
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

cs = ax1.contourf(x, y, wn0, levels=60, cmap="YlGnBu")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_title("t=0")
#fig.colorbar(cs, ax = ax1)

cs = ax2.contourf(x, y, wn10, levels=60, cmap="YlGnBu")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_title("t=10")
#fig.colorbar(cs, ax = ax2)

cs = ax3.contourf(x, y, wn20, levels=60, cmap="YlGnBu")
ax3.set_xlabel("X")
ax3.set_ylabel("Y")
ax3.set_title("t=20")
#fig.colorbar(cs, ax = ax3)

fig.tight_layout()

fig.subplots_adjust(bottom=0.2)
cbaxes = fig.add_axes([0.2, 0.05, 0.6, 0.05])
fig.colorbar(cs, orientation="horizontal", cax = cbaxes)

filename = "vortex-merger"+str(nx)+".pdf"
fig.savefig( filename )
plt.show()

