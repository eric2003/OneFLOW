import numpy as np
import matplotlib.pyplot as plt
    
def plot_all_cell_center( xcc, yref ):
    plt.scatter(xcc, np.full_like(xcc, yref), s=20, facecolor='black', edgecolor='black', linewidth=1)
    return
    
def plot_cell_center( xcc, yref ):
    nx = xcc.size
    ii = nx // 2
    im = ii - 1
    ip = ii + 1
    xcc_new = []
    for i in range(0, nx):
        if i > 1 and i < im:
            continue
        if i > ip and i < nx-2:
            continue
        xcc_new.append( xcc[i] )
    plt.scatter(xcc_new, np.full_like(xcc_new, yref), s=20, facecolor='black', edgecolor='black', linewidth=1)
    return
    
def plot_ghost_cell_center( xcc, yref ):
    plt.scatter(xcc, np.full_like(xcc, yref), s=20, facecolor='red', edgecolor='black', linewidth=1)
    return
    
def plot_ghost_mesh( x, yref ):
    dx = x[1] - x[0]
    dy = 0.1 * dx
    for xm in x:
        plt.plot([xm, xm], [yref-dy, yref+dy], 'k-')  # 绘制垂直线
    return

def plot_mesh( x, yref ):
    dx = x[1] - x[0]
    dy = 0.1 * dx
    for xm in x:
        plt.plot([xm, xm], [yref-dy, yref+dy], 'k-')  # 绘制垂直线
    #
    
    nxc = x.size - 1
    ii = nxc // 2
    im = ii - 1
    ip = ii + 1
    
    for i in range(0, nxc):
        if i > 1 and i < im:
            plt.plot([x[i], x[i+1]], [yref, yref], 'k--', linewidth=1)
        elif i > ip and i < nx-2:
            plt.plot([x[i], x[i+1]], [yref, yref], 'k--', linewidth=1)
        else :
            plt.plot([x[i], x[i+1]], [yref, yref], 'b-', linewidth=1)
    #plt.plot(x, np.full_like(x, yref), 'k--', linewidth=1)
    return
    
def genstrNp(i):
    if i != 0:
        ai = abs(i)
        if i > 0:
           ss = '+'
        else:
           ss = '-'
        mystr = r'$N' + ss + f'{abs(i)}'+r'$'
    else:
        mystr = r'$N'+r'$'
    return mystr
    
def genstrXNhalf(strn, i):
    if i != 0:
        ai = abs(i)
        if i > 0:
           ss = '+'
        else:
           ss = '-'
        mystr = r'$x_{' + strn + ss + r'\frac{' + f'{ai}' + r'}{2}}$'
    else:
        mystr = r'$x_{' + strn + ss + r'\frac{' + r'}{2}}$'
    return mystr    
   
def plot_label(x, xcc, yref, ishift):
    x0 = x[0]
    dx = x[1] - x[0]
    dyb = 0.8 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    ybc = yref - 0.5* dyb
    
    i1 = 2 * ( 1 + ishift ) - 1
    i2 = i1 + 2
    
    str1 = genstrXNhalf('', i1)
    str2 = genstrXNhalf('', i2)
    
    plt.text(x[0], yb, str1, fontsize=12, ha='center')
    plt.text(x[1], yb, str2, fontsize=12, ha='center')
    
    i1 = 2 * ( 1 + ishift ) - 1
    i2 = i1 - 2
    
    str1 = genstrXNhalf('N', i1)
    str2 = genstrXNhalf('N', i2)
    
    plt.text(x[-1], yb, str1, fontsize=12, ha='center')
    plt.text(x[-2], yb, str2, fontsize=12, ha='center')
    
    
    plt.text(x[0], yt, r'$x=a$', fontsize=12, ha='center')
    plt.text(x[-1], yt, r'$x=b$', fontsize=12, ha='center')
    
    i1 = 1 + ishift
    i2 = i1 + 1
    
    str1 = r'$'+f'{i1}'+r'$'
    str2 = r'$'+f'{i2}'+r'$'
    
    plt.text(xcc[0], ybc, str1, fontsize=12, ha='center')
    plt.text(xcc[1], ybc, str2, fontsize=12, ha='center')
    
    i1 = ishift
    i2 = i1 - 1
    
    str1 = genstrNp(i1)
    str2 = genstrNp(i2)
    
    plt.text(xcc[-1], ybc, str1, fontsize=12, ha='center')
    plt.text(xcc[-2], ybc, str2, fontsize=12, ha='center')    
    
    nx = xcc.size
    i = nx // 2
    print("i=",i)
    im = i - 1
    ip = i + 1
    
    plt.text(xcc[im], ybc, r'$i-1$', fontsize=12, ha='center')
    plt.text(xcc[i], ybc, r'$i$', fontsize=12, ha='center')
    plt.text(xcc[ip], ybc, r'$i+1$', fontsize=12, ha='center')
    
    str = r'$a=x_{-\frac{1}{2}}<x_{\frac{1}{2}}<\cdots<x_{N-\frac{3}{2}}<x_{N-\frac{1}{2}}=b$'
    str = 'Grid: ' + str
    
    nx = xcc.size
    ii = nx // 2    
    
    plt.text(x[ii], yb-dx, str, fontsize=12, ha='center')
 
    return
    
def plot_ghost_label_left(xg, xgcc, yref, ishift):
    dx = abs(xg[1] - xg[0])
    dyb = 0.8 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    ybc = yref - 0.5* dyb
    
    i1 = -1 + 2 * ishift
    i2 = i1 - 2
    
    str1 = r'$x_{-\frac{' + f'{abs(i1)}' + r'}{2}}$'
    str2 = r'$x_{-\frac{' + f'{abs(i2)}' + r'}{2}}$'
    
    plt.text(xg[1], yb, str1, fontsize=12, ha='center')
    plt.text(xg[2], yb, str2, fontsize=12, ha='center')
    
    i1 = 0+ishift
    i2 = i1-1
    str1 = r'$'+f'{i1}'+r'$'
    str2 = r'$'+f'{i2}'+r'$'
    
    plt.text(xgcc[0], ybc, str1, fontsize=12, ha='center')
    plt.text(xgcc[1], ybc, str2, fontsize=12, ha='center')
 
    return

   
def plot_ghost_label_right(xg, xgcc, yref, ishift):
    dx = abs(xg[1] - xg[0])
    dyb = 0.8 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    ybc = yref - 0.5* dyb
    
    i1 = 2*(1+ishift)+1
    i2 = i1 + 2
    
    str1 = r'$x_{N+\frac{' + f'{abs(i1)}' + r'}{2}}$'
    str2 = r'$x_{N+\frac{' + f'{abs(i2)}' + r'}{2}}$'
    plt.text(xg[1], yb, str1, fontsize=12, ha='center')
    plt.text(xg[2], yb, str2, fontsize=12, ha='center')
    
    i1 = 1+ishift
    i2 = i1+1
    str1 = genstrNp(i1)
    str2 = genstrNp(i2)
    
    plt.text(xgcc[0], ybc, str1, fontsize=12, ha='center')
    plt.text(xgcc[1], ybc, str2, fontsize=12, ha='center')    
 
    return

# 设置字体为 Times New Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

# 设置图形大小和样式
plt.figure(figsize=(12, 5))

nx = 9
L  = 1.0
x_l = 0.0
dx = L / nx

x   = np.zeros(nx+1, dtype=np.float64)
xcc = np.zeros(nx, dtype=np.float64)

nghost = 2
x_ghost_l   = np.zeros(nghost+1, dtype=np.float64)
xcc_ghost_l = np.zeros(nghost, dtype=np.float64)
x_ghost_r   = np.zeros(nghost+1, dtype=np.float64)
xcc_ghost_r = np.zeros(nghost, dtype=np.float64)

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)
    
for i in range(0, nx):
    xcc[i] = 0.5*(x[i]+x[i+1])
    
x_ghost_l[0] = x[0]
for ighost in range(1, nghost+1):
    dx = x[0] - x[ighost]
    x_ghost_l[ighost] = x[0] + dx
    
x_ghost_r[0] = x[nx]
for ighost in range(1, nghost+1):
    dx = x[nx] - x[nx-ighost]
    x_ghost_r[ighost] = x[nx] + dx
    
for ighost in range(0, nghost):
    xcc_ghost_l[ighost] = 0.5*(x_ghost_l[ighost]+x_ghost_l[ighost+1])
    xcc_ghost_r[ighost] = 0.5*(x_ghost_r[ighost]+x_ghost_r[ighost+1])
    
print("x=",x)
print("xcc=",xcc)

print("x_ghost_l=",x_ghost_l)
print("xcc_ghost_l=",xcc_ghost_l)
print("x_ghost_r=",x_ghost_r)
print("xcc_ghost_r=",xcc_ghost_r)

yref = 0.0

plot_ghost_cell_center( xcc_ghost_l,yref )
plot_ghost_cell_center( xcc_ghost_r,yref )

plot_ghost_mesh( x_ghost_l, yref )
plot_ghost_mesh( x_ghost_r, yref )

ishift = -1
plot_ghost_label_left(x_ghost_l, xcc_ghost_l, yref, ishift)
plot_ghost_label_right(x_ghost_r, xcc_ghost_r, yref, ishift)

plot_cell_center( xcc, yref )
plot_mesh( x, yref )
plot_label(x, xcc, yref, ishift)
    
plt.axis('equal')
plt.axis('off')

plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
plt.show()