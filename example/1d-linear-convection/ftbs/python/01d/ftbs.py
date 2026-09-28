import numpy as np
import matplotlib.pyplot as plt
import time, sys                   #and load some utilities

def set_fuction(u, ct) :
    for i in range(nx):
        xm = i * dx
        if 0.5 <= xm - ct <= 1:
            u[i] = 2
        else:
            u[i] = 1

nx = 41  # try changing this number from 41 to 81 and Run All ... what happens?
dx = 2 / (nx-1)
nt = 25    #nt is the number of timesteps we want to calculate
dt = .025  #dt is the amount of time each timestep covers (delta t)
c = 1      #assume wavespeed of c = 1

total_t = dt * nt
utheory = np.ones(nx)
x = np.linspace(0, 2, nx)
u = np.ones(nx)

set_fuction(u,0)
set_fuction(utheory,c * total_t)
       
un = np.ones(nx) #initialize a temporary array

for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
    un = u.copy() ##copy the existing values of u into un
    for i in range(1, nx): ## you can try commenting this line and...
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
        
plt.plot(x, utheory, "k-", linewidth=1.0, label="Exact solution")
plt.scatter(x, u, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="FTBS scheme")
plt.xlabel("$x$")
plt.ylabel("$u$")
plt.title("Solution field")
plt.legend()
plt.tight_layout()
plt.show();