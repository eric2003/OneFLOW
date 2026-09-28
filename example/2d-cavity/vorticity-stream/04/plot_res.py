import numpy as np
import matplotlib.pyplot as plt

iter_hist = []
res_hist = []

with open('residual_plot.txt', 'r') as f:
    for index, line in enumerate(f):
        words = line.strip().split()
        iter_hist.append( float(words[0]) )
        res_hist.append( float(words[1]) )

print("iter_hist=\n",iter_hist)
print("res_hist=\n",res_hist)
        
fig = plt.figure("OneFLOW CFD Solver", figsize=(8, 6), dpi=100)

ax = fig.add_subplot(1,1,1)
ax.set_yscale("log", base=10)
#ax.set_xlim(0,10000)
ax.set_xlabel("Iteration count")
ax.set_ylabel(r"Residual $|L|_2$ norm")
ax.plot(iter_hist, res_hist, color="red", linewidth=2.0, label="rms")
ax.legend() 
fig.tight_layout()
fig.savefig("ldc_residual.pdf")
plt.show()


