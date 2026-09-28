import numpy as np
import matplotlib.pyplot as plt

n = 12
m = int(n/2)
dx = 1.0 / n
h = dx

xh = np.zeros( n + 1 )
yh = np.zeros( n + 1 )
ih = np.zeros( n + 1 )

xh_odd = np.zeros( m + 1 )
yh_odd = np.zeros( m + 1 )
ih_odd = np.zeros( m + 1 )

xh_even = np.zeros( m + 1 )
yh_even = np.zeros( m + 1 )
ih_even = np.zeros( m + 1 )

x2h = np.zeros( m + 1 )
y2h = np.zeros( m + 1 )
i2h = np.zeros( m + 1 )

km = 2
io = 0
ie = 0
for j in range(0, n + 1):
    xh[j] = j * h
    ih[j] = j
    yh[j] = np.sin( j * km * np.pi / n )
    if (j % 2) == 0:
        xh_odd[io] = xh[j]
        ih_odd[io] = ih[j]
        yh_odd[io] = yh[j]
        io += 1
    else:
        xh_even[ie] = xh[j]
        ih_even[ie] = ih[j]
        yh_even[ie] = yh[j]
        ie += 1

print("ih_odd=\n",ih_odd)
print("ih_even=\n",ih_even)

for j in range(0, m + 1):
    i2h[j] = j
    x2h[j] = j * 2 * h
    y2h[j] = np.sin( j * km * np.pi / m )
fig = plt.figure('OneFLOW-CFD Solver Wave with Wavenumber', figsize=(6, 4), dpi=100)

ax = fig.add_subplot(2,1,1)
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel("X")
ax.set_ylabel("Y")
mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(i2h,y2h,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)
for j in range(0, m+1):
    ax.vlines(x = i2h[j], ymin = -1.2, ymax = 1.2, linestyle = ':')
#ax.legend()

ax = fig.add_subplot(2,1,2)
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel("X")
ax.set_ylabel("Y")
mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'fine grid'
ax.plot(ih,yh,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)
ax.plot(ih_odd,yh_odd,color=mycolor,linewidth=1,linestyle = 'none',marker=mymarker,markerfacecolor='none',label=mylabel)
ax.plot(ih_even[:ie],yh_even[:ie],color=mycolor,linewidth=1,linestyle = 'none',marker=mymarker,label=mylabel)
ax.plot(ih,yh,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)

for j in range(0, n+1):
    ax.vlines(x = ih[j], ymin = -1.2, ymax = 1.2, linestyle = ':')
#ax.legend()

fig.tight_layout()
plt.show()
