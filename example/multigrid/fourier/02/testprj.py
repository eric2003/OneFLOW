import numpy as np
import matplotlib.pyplot as plt
    
n = 12
h = 1.0/n
dx = h

print("dx=",dx)

x = np.zeros( n+1 )
kk=np.array([1,2,3,4,6,8,9])
num = len(kk)
print("num=",num)
yy = np.zeros( (num,n+1) )

for i in range(0, n+1):
    x[i] = i * dx

for j in range(0, num):
  for i in range(0, n+1):
      yy[j][i] = np.sin( i * kk[j] * np.pi / n )
    
id = np.zeros(n+1)
for i in range(0, n+1):
   id[i] = i
   
fig = plt.figure("OneFLOW-CFD Solver Graphs of the Fourier modes", figsize=(6, 10), dpi=100)  

for j in range(0, num):
    ax = fig.add_subplot(num,1,j+1)
    mylabel="k = {value}".format( value = kk[j])
    ax.plot(id,yy[j],color='black',linewidth=1,marker='o',markerfacecolor='none',label=mylabel)
    ax.legend()

fig.tight_layout()
plt.show()
