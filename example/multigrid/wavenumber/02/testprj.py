import numpy as np
import matplotlib.pyplot as plt

n = 12
m = int(n/2)
dx = 1.0/n
h = dx
nPoints=np.array([n,n/2])
mylabel_list=np.array(["k=4 wave on n = 12 grid","k=4 wave on n = 6 grid"])
mycolor_list=np.array(['blue','red'])
mymarker_list=np.array(['o','s'])

xh = np.zeros( n + 1 )
yh = np.zeros( n + 1 )

x2h = np.zeros( m + 1 )
y2h = np.zeros( m + 1 )

km = 4

nf = 360
xf = np.zeros( nf + 1 )
yf = np.zeros( nf + 1 )
dxf= 1.0/nf

for j in range(0, nf + 1):
    xf[j] = j * dxf
    yf[j] = np.sin( j * km * np.pi / nf )

for j in range(0, n + 1):
    xh[j] = j * h
    yh[j] = np.sin( j * km * np.pi / n )


for j in range(0, m + 1):
    x2h[j] = j * 2 * h
    y2h[j] = np.sin( j * km * np.pi / m )
fig = plt.figure('OneFLOW-CFD Solver Wave with Wavenumber', figsize=(6, 8), dpi=100)

ax = fig.add_subplot(3,1,1)
ax.set_xlabel("X")
ax.set_ylabel("Y")
mycolor = mycolor_list[0]
mymarker = mymarker_list[0]
mylabel = mylabel_list[0]
ax.plot(xf,yf,color='black',linewidth=1, marker='none')
ax.plot(xh,yh,color=mycolor,linewidth=1,linestyle='none', marker=mymarker,markerfacecolor='none',label=mylabel)
mycolor = mycolor_list[1]
mymarker = mymarker_list[1]
mylabel = mylabel_list[1]
ax.plot(x2h,y2h,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)
ax.legend()

ax = fig.add_subplot(3,1,2)
ax.set_xlabel("X")
ax.set_ylabel("Y")
mycolor = mycolor_list[0]
mymarker = mymarker_list[0]
mylabel = mylabel_list[0]
ax.plot(xh,yh,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)
ax.legend()

ax = fig.add_subplot(3,1,3)
ax.set_xlabel("X")
ax.set_ylabel("Y")
mycolor = mycolor_list[1]
mymarker = mymarker_list[1]
mylabel = mylabel_list[1]
ax.plot(x2h,y2h,color=mycolor,linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)
ax.legend()

fig.tight_layout()
plt.show()
