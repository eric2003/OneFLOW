# LIBRARY
# vector manipulation
import numpy as np
# math functions
import math 

# THIS IS FOR PLOTTING
import matplotlib.pyplot as plt # side-stepping mpl backend
import warnings
warnings.filterwarnings("ignore")

N=20
Nt=10
h=2*np.pi/N
k=1/Nt
r=k/(h*h)
time_steps=10
time=np.arange(0,(time_steps+.5)*k,k)
x=np.arange(0,2*np.pi+h/2,h)

w=np.zeros((N+1,time_steps+1))
b=np.zeros(N-1)
# Initial Condition
for i in range (0,N+1):
    w[i,0]=1-np.cos(x[i])
    
ipos = np.zeros(N+1)
ineg = np.zeros(N+1)
for i in range(0,N+1):    
   ineg[i] = i-1
   ipos[i] = i+1
ipos[N] = 0
ineg[0] = N

lamba=k/h
for j in range(0,time_steps):
    for i in range (0,N+1):
        w[i,j+1]=(w[int(ipos[i]),j]+w[int(ineg[i]),j])/2-lamba/2*(w[int(ipos[i]),j]-w[int(ineg[i]),j])

fig = plt.figure(figsize=(12,6))
plt.plot(x,w[:,Nt-1],'o:')
plt.xlabel('x')
plt.ylabel('w')

fig.tight_layout()
plt.show()