import numpy as np
import matplotlib.pyplot as plt

def Gauss_Seidel(v, niter, dx):
    dx2 = dx * dx
    den = - 2.0 / dx2
    res = np.zeros(niter)
    for iter in range(0, niter):
        for i in range(1, n):
            d2udx2 = ( v[i-1] - 2*v[i] + v[i+1] )/dx2
            rm = -d2udx2
            v[i] = v[i] + rm / den
        resm = np.abs(np.max(v))
        res[iter] = resm
    return res
    
n = 12
h = 1.0/n
dx = h

print("dx=",dx)

x = np.zeros( n+1 )
y1 = np.zeros( n+1 )
y2 = np.zeros( n+1 )
y3 = np.zeros( n+1 )
y4 = np.zeros( n+1 )

k1=1
k2=2
k3=3
k4=4

kk=np.array([1,2,3,4,6,8,9])
num = len(kk)
yy = np.zeros( (num,n+1) )

for i in range(0, n+1):
    x[i] = i * dx
    y1[i] = np.sin( i * k1 * np.pi / n )
    y2[i] = np.sin( i * k2 * np.pi / n )
    y3[i] = np.sin( i * k3 * np.pi / n )
    y4[i] = np.sin( i * k4 * np.pi / n )
    

id = np.zeros(n+1)
for i in range(0, n+1):
   id[i] = i
   
#print("v=\n",v)   

fig = plt.figure("OneFLOW-CFD Solver Graphs of the Fourier modes", figsize=(6, 4), dpi=100)  

ax1 = fig.add_subplot(411)
ax1.plot(id,y1,color='black',linewidth=1,marker='o',markerfacecolor='none',label='k=1')
#ax1.legend()

ax2 = fig.add_subplot(412)
ax2.plot(id,y2,color='black',linewidth=1,marker='o',markerfacecolor='none',label='k=2')
#ax2.legend()

ax3 = fig.add_subplot(413)
ax3.plot(id,y3,color='black',linewidth=1,marker='o',markerfacecolor='none',label='k=3')
#ax3.legend()

ax4 = fig.add_subplot(414)
ax4.plot(id,y4,color='black',linewidth=1,marker='o',markerfacecolor='none',label='k=4')
#ax4.legend()
fig.tight_layout()
plt.show()
