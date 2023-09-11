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

plt.axis('off')

txt1 = "}"
plt.text(x[0]+0.1, y[0]+0.1, txt1, fontsize=50)

plt.text(x[1], y[1]-0.1, txt1, fontsize=100, rotation=-90, transform_rotates_text=True)

for j in range(n+1):
    plt.hlines(y = y[j], xmin = x[0], xmax = x[-1], colors="k" )
    for i in range(n+1):
        text = "({0},{1})".format(i,j)
        plt.text(x[i], y[j], text, fontsize=10)
        plot_point( x[i], y[j] )

for i in range(n+1):
    plt.vlines(x = x[i], ymin = y[0], ymax = y[-1], colors="k")
    
plt.show()