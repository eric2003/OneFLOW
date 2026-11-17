import numpy as np
import matplotlib.pyplot as plt
    
def plot_cell_center( xcc, yref ):
    plt.scatter(xcc, np.full_like(xcc, yref), s=10, facecolor='black', edgecolor='black', linewidth=1)
    return

def plot_mesh( x, yref ):
    dx = x[1] - x[0]
    dy = 0.1 * dx
    for xm in x:
        plt.plot([xm, xm], [yref-dy, yref+dy], 'k-')  # 绘制垂直线
    #
    plt.plot(x, np.full_like(x, yref), 'k--', linewidth=1)
    return
    
   
def plot_label(x, xcc, yref):
    x0 = x[0]
    dx = x[1] - x[0]
    dyb = 0.5 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    plt.text(x[0], yb, r'$x_{\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[1], yb, r'$x_{\frac{3}{2}}$', fontsize=12, ha='center')
    plt.text(x[-2], yb, r'$x_{N-\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[-1], yb, r'$x_{N+\frac{1}{2}}$', fontsize=12, ha='center')
    
    plt.text(x[0], yt, r'$x=0$', fontsize=12, ha='center')
    plt.text(x[-1], yt, r'$x=L$', fontsize=12, ha='center')
    
    return    

# 设置字体为 Times New Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

# 设置图形大小和样式
plt.figure(figsize=(12, 5))

nx = 10
L  = 1.0
x_l = 0.0
dx = L / nx

x   = np.zeros(nx+1, dtype=np.float64)
xcc = np.zeros(nx, dtype=np.float64)

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)
    
for i in range(0, nx):
    xcc[i] = 0.5*(x[i]+x[i+1])
    
print("x=",x)
print("xcc=",xcc)

yref = 0.0

plot_cell_center( xcc, yref )
plot_mesh( x, yref )
plot_label(x, xcc, yref)
    
plt.axis('equal')
plt.axis('off')

plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
plt.show()