#https://ikuz.eu/machine-learning-and-computer-science/the-concept-of-conjugate-gradient-descent-in-python/

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

A = np.matrix([[3.0, 2.0], [2.0, 6.0]])
b = np.matrix([[2.0], [-8.0]])  # we will use the convention that a vector is a column vector
c = 0.0

def f(x, A, b, c):
    return np.sum(0.5 * x.T * A * x - b.T * x + c)
    
def bowl(A, b, c):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(projection='3d')
    size = 20
    x1 = list(np.linspace(-6, 6, size))
    x2 = list(np.linspace(-6, 6, size))
    x1, x2 = np.meshgrid(x1, x2)
    zs = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            x = np.matrix([[x1[i,j]], [x2[i,j]]])
            zs[i,j] = f(x, A, b, c)
    ax.plot_surface(x1, x2, zs, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0)
    plt.show()
    return x1, x2, zs    
    

x1, x2, zs = bowl(A, b, c)    

def contoursteps(x1, x2, zs, steps=None):
    fig = plt.figure(figsize=(6,6))
    cp = plt.contour(x1, x2, zs, 10)
    plt.clabel(cp, inline=1, fontsize=10)
    if steps is not None:
        steps = np.matrix(steps)
        plt.plot(steps[:,0], steps[:,1], '-o')
    plt.show()
    
#The Method of Steepest Descent    
x = np.matrix([[-2.0],[-2.0]])
steps = [(-2.0, -2.0)]

i = 0
imax = 10
eps = 0.01
r = b - A * x
delta = np.sum(r.T * r)
delta0 = delta
    
while i < imax and delta > eps**2 * delta0:
    alpha = delta / np.sum(r.T * (A * r))
    x = x + alpha * r
    steps.append((x[0,0], x[1,0]))  # store steps for future drawing
    r = b - A * x
    delta = np.sum(r.T * r)
    i += 1    
print("niter=",i)    
    
contoursteps(x1, x2, zs, steps)    
