import numpy as np
import matplotlib.pyplot as plt

iter_hist = []
res_hist = []

with open('cg_residual.txt', 'r') as f:
    for index, line in enumerate(f):
        words = line.strip().split()
        iter_hist.append( float(words[0]) )
        res_hist.append( float(words[2]) )
        
fig = plt.figure("OneFLOW CFD Poisson Solver", figsize=(8, 6), dpi=100)

ax = fig.add_subplot(1,1,1)
ax.set_yscale("log", base=10)
ax.set_xlabel("Iteration count")
ax.set_ylabel(r"$|r|_2$")
ax.plot(iter_hist, res_hist, color="blue", linewidth=4.0, label="Conjugate-Gradient method")

ax.legend() 
fig.tight_layout()
plt.show()
