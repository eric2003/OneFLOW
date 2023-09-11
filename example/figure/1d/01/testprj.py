import matplotlib.pyplot as plt
import numpy as np

def plot_point( xp, yp ):
    x = [xp]
    y = [yp]
    plt.plot(x, y, marker="o", markersize=8, markeredgecolor="blue",
    markerfacecolor="blue")

n=10
h=1/n
print("n=",n)
print("h=",h)

x = np.zeros(n+1)
xmin = 0
xmax = 1
ymin = -0.02
ymax = 0.02
x[0] = xmin
for i in range(n):
    x[i+1] = xmin + ( i + 1 ) * h
    

plt.axis('off')

plt.ylim([-1, 1])

plt.hlines(y = 0, xmin = x[0], xmax = x[-1], colors="k" )
for i in range(n+1):
    text = "({0},{1})".format(i,0)
    #plt.text(x[i], 0, text, fontsize=10)
    #plot_point( x[i], 0)
tleft=r'$\Omega^{h}:$'
tb0=r'$x=0$'
tb1=r'$x=1$'
tx0=r'$x_{0}$'
tx1=r'$x_{1}$'
tx2=r'$x_{2}$'
txj=r'$x_{j}$'
txn1=r'$x_{n-1}$'
txn=r'$x_{n}$'
for i in range(n+1):
    plt.vlines(x = x[i], ymin = ymin, ymax = ymax, colors="k")
    

plt.text(x[0]-1.2*h, -0.05, tleft, fontsize=15)    
plt.text(x[0]-0.25*h, 0.1, tb0, fontsize=10)
plt.text(x[-1]-0.25*h, 0.1, tb1, fontsize=10)
plt.text(x[0]-0.25*h, -0.1, tx0, fontsize=10)
plt.text(x[1]-0.25*h, -0.1, tx1, fontsize=10)
plt.text(x[2]-0.25*h, -0.1, tx2, fontsize=10)
plt.text(x[5]-0.25*h, -0.1, txj, fontsize=10)
plt.text(x[-2]-0.25*h, -0.1, txn1, fontsize=10)
plt.text(x[-1]-0.25*h, -0.1, txn, fontsize=10)
    
plt.show()