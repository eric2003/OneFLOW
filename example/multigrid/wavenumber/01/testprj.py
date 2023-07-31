import numpy as np
import matplotlib.pyplot as plt

n = 361
x = np.zeros( n )
y1 = np.zeros( n )
y3 = np.zeros( n )
y6 = np.zeros( n )

for i in range(0, n):
    x[i] = i
    y1[i] = np.sin( i * 1 * np.pi / n )
    y3[i] = np.sin( i * 3 * np.pi / n )
    y6[i] = np.sin( i * 6 * np.pi / n )


plt.figure("OneFLOW-CFD Solver Wavenumber", figsize=(6, 4), dpi=100)  
plt.plot(x,y1,'k-',linewidth=2,label='k=1')
plt.plot(x,y3,'b-.',linewidth=2,label='k=3')
plt.plot(x,y6,'g--',linewidth=2,label='k=6')
plt.legend()
plt.show()
