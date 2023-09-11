import matplotlib.pyplot as plt
import numpy as np

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

for j in range(n+1):
    plt.hlines(y = y[j], xmin = x[0], xmax = x[-1])

for i in range(n+1):
    plt.vlines(x = x[i], ymin = y[0], ymax = y[-1])
    
plt.show()