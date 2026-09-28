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
        if i > 0 and i < im:
            continue
        if i > ip and i < nx-1:
            continue
        xcc_new.append( xcc[i] )
    plt.scatter(xcc_new, np.full_like(xcc_new, yref), s=20, facecolor='black', edgecolor='black', linewidth=1)
    return 

def plot_cell_center_rs( xcc, yref, r, s ):
    nx = xcc.size
    ii = nx // 2
    xcc_new = []
    for m in range(-r, s+1):
        xcc_new.append( xcc[ii+m] )
    plt.scatter(xcc_new, np.full_like(xcc_new, yref), s=20, facecolor='black', edgecolor='black', linewidth=1)
    return    

def plot_mesh( x, yref ):
    dx = x[1] - x[0]
    dy = 0.1 * dx
    nx = x.size
    for i in range(0, nx):
        xm = x[i]
        plt.plot([xm, xm], [yref-dy, yref+dy], 'k-')  # 绘制垂直线
    
    nxc = x.size - 1
    
    for i in range(0, nxc):
        plt.plot([x[i], x[i+1]], [yref, yref], 'b-', linewidth=1)
    return
    
def plot_mesh_rs( x, yref, r, s):
    dx = x[1] - x[0]
    dy = 0.1 * dx
   
    nxc = xcc.size
    ii = nxc // 2
    
    idc = []
    
    for m in range(-r, s+1):
        #print(f'm={m}, r={r}, s={s}')
        idc.append( ii+m )
    #print(f"idc={idc}")
    idv = idc.copy()
    idv.append(idc[-1]+1)
    
    ncell = len( idc )
    nvertex = len( idv )
    for i in range(0, nvertex):
        xm = x[ idv[i] ]
        plt.plot([xm, xm], [yref-dy, yref+dy], 'k-')  # 绘制垂直线
        
    for i in range(0, ncell):
        plt.plot([x[idc[i]], x[idc[i]+1]], [yref, yref], 'b-', linewidth=1)
    return    
    
def plot_label(x, xcc, yref):
    dx = x[1] - x[0]
    dyb = 0.5 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    ybc = yref - 0.5* dyb
    plt.text(x[0], yb, r'$x_{i-\frac{5}{2}}$', fontsize=12, ha='center')
    plt.text(x[1], yb, r'$x_{i-\frac{3}{2}}$', fontsize=12, ha='center')
    plt.text(x[2], yb, r'$x_{i-\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[3], yb, r'$x_{i+\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[4], yb, r'$x_{i+\frac{3}{2}}$', fontsize=12, ha='center')
    plt.text(x[5], yb, r'$x_{i+\frac{5}{2}}$', fontsize=12, ha='center')
    
    nx = xcc.size
    i = nx // 2
    print("i=",i)
    im  = i - 1
    im1 = i - 2
    ip  = i + 1
    ip1 = i + 2
    
    plt.text(xcc[im1], ybc, r'$i-2$', fontsize=12, ha='center')
    plt.text(xcc[im], ybc, r'$i-1$', fontsize=12, ha='center')
    plt.text(xcc[i], ybc, r'$i$', fontsize=12, ha='center')
    plt.text(xcc[ip], ybc, r'$i+1$', fontsize=12, ha='center')
    plt.text(xcc[ip1], ybc, r'$i+2$', fontsize=12, ha='center')
    return
    
def plot_label_rs(x, xcc, yref, r, s):
    dx = x[1] - x[0]
    dyb = 0.5 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    ybc = yref - 0.5* dyb
    ytt = yref + 0.5* dyt
    namelist = []
    namelist.append(r'$x_{i-\frac{5}{2}}$')
    namelist.append(r'$x_{i-\frac{3}{2}}$')
    namelist.append(r'$x_{i-\frac{1}{2}}$')
    namelist.append(r'$x_{i+\frac{1}{2}}$')
    namelist.append(r'$x_{i+\frac{3}{2}}$')
    namelist.append(r'$x_{i+\frac{5}{2}}$')
    
    nx = xcc.size
    ii = nx // 2
    
    idc = []
    for m in range(-r, s+1):
        idc.append( ii+m )
    idv = idc.copy()
    idv.append(idc[-1]+1)
    
    ncell = len( idc )
    nvertex = len( idv )
    
    for i in range(0, nvertex):
        xm = x[ idv[i] ]
        name = namelist[ idv[i] ]
        #plt.text(xm, yb, name, fontsize=12, ha='center')
        plt.text(xm, ytt, name, fontsize=12, ha='center')
    
    for m in range(-r, s+1):
        ss = '-'
        if m > 0 :
           ss = '+'
        str = r'$i' + ss + f'{abs(m)}' + r'$'
        if m == 0 :
            plt.text(xcc[ii+m], ybc, r'$i$', fontsize=12, ha='center')
        else:
            plt.text(xcc[ii+m], ybc, str, fontsize=12, ha='center')
            
    str = r'$' + f'({r=},{s=})' + r'$'
    ishift = (-r+s)//2
    plt.text(xcc[ii+ishift], yb, str, fontsize=12, ha='center')
    
    return
    
def getrs(k,rv,sv):
    kk = k-1
    for m in range(0, k):
        s = m
        r = kk - s
        rv.append( r )
        sv.append( s )
    return    

# 设置字体为 Times New Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

# 设置图形大小和样式
plt.figure(figsize=(12, 6))

nx = 5
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

#k=3
k=2
rv = []
sv = []
getrs(k,rv,sv)
print(f'{rv=},{sv=}')

dyref = 0.2

size = len(rv)
print(f'{size=}')
for i in range(0, size):
    yref = 0.0 - i * dyref
    r=rv[i]
    s=sv[i]
    plot_cell_center_rs( xcc, yref, r, s)
    plot_mesh_rs( x, yref, r, s)
    plot_label_rs(x, xcc, yref, r, s)

    
plt.axis('equal')
plt.axis('off')

plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
plt.show()