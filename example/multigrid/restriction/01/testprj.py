import matplotlib.pyplot as plt
import numpy as np

def plot_point( xp, yp ):
    x = [xp]
    y = [yp]
    plt.plot(x, y, marker="o", markersize=8, markeredgecolor="blue",
    markerfacecolor="blue")

n=4
h=1/n
print("n=",n)
print("h=",h)

x = np.zeros(n+1)
y = np.zeros(n+1)
xmin = 0
xmax = 1
ymin = 0
ymax = 1
x[0] = xmin
for i in range(n):
    #print("i=",i,"n=",n)
    x[i+1] = xmin + ( i + 1 ) * h

for j in range(n):
    y[j+1] = ymin + ( j + 1 ) * h
    
ue = np.zeros( ( n + 1, n + 1 ) )

for j in range(1, n):
    for i in range(0, n+1):
        y16m = np.sin( i * 16 * np.pi / n )
        y40m = np.sin( i * 40 * np.pi / n )
        ue[i,j] = 0.5 * ( y16m + y40m )
    
plt.figure("OneFLOW-CFD Solver+Restriction", figsize=(6, 4), dpi=100)
plt.axis('off')
#plot_point( x[0], y[0] )

plt.subplot(1, 2, 1)

plt.title( "Line Value" )

for j in range(n+1):
    print("j=",j,"ue=\n",ue[:,j])
    plt.plot(x,ue[:,j],color='black',linewidth=1)

plt.subplot(1, 2, 2)

for j in range(n+1):
    plt.hlines(y = y[j], xmin = x[0], xmax = x[-1], colors="k" )
    for i in range(n+1):
        text = "({0},{1})".format(i,j)
        plt.text(x[i], y[j], text, fontsize=10)
        plot_point( x[i], y[j] )

for i in range(n+1):
    plt.vlines(x = x[i], ymin = y[0], ymax = y[-1], colors="k")
    

plt.contourf(x, y, ue, levels=20, cmap="jet")
#plt.contour(x, y, ue, levels=20, cmap="jet")
#plt.colorbar()
plt.title( "Exact solution" )
plt.tight_layout()
plt.show()
