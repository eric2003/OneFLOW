import numpy as np
import matplotlib.pyplot as plt
    
n = 48
h = 1.0/n
dx = h

print("dx=",dx)

x = np.zeros( n+1 )
eigen = np.zeros( n+1 )

for i in range(0, n+1):
    x[i] = i * dx

for k in range(0, n+1):
  value = np.cos( k * np.pi / n )
  v2 = value * value
  eigen[k] = v2
    
id = np.zeros(2*n+1)
for i in range(0, 2*n+1):
   id[i] = i
   
fig = plt.figure("OneFLOW-CFD Solver Eigenvalues of Gauss–Seidel iteration matrix", figsize=(6, 4), dpi=100)  
ax = fig.add_subplot()

#ax.spines['left'].set_position('center')
#ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel("k")
ax.set_ylabel(r'$\lambda_{k}(P_{G})$',rotation='horizontal')
ax.yaxis.set_label_coords(0.0, 1.0)
ax.xaxis.set_label_coords(1.01, 0.0)

ax.plot(x,eigen,color='blue',linewidth=1,marker='o',markerfacecolor='none')

fig.tight_layout()
plt.show()
