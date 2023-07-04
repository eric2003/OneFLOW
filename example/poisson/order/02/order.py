import numpy as np
import matplotlib.pyplot as plt

grid = [32,64,128,256,512]
fft_spectral = [1.339154672220572e-16, 1.342489800590083e-16, 1.3269947532677886e-16, \
                1.4514549677756118e-16 , 1.485803753784319e-16]
fft_fdm = [0.0015607100315532957 , 0.0005987381110678801 , 0.00014313734718665358, \
           3.549617203207291e-5, 8.865373334924762e-6]
ref_x = [64,256]
ref_y = [2e-3, 1.25e-4]

fig = plt.figure("OneFLOW-CFD Solver+Order", figsize=(10, 4), dpi=100)

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.plot(grid, fft_spectral, color="red",lw=4,
                  marker = "o", markeredgecolor="k",
                  markersize=12)
ax1.set_xscale("log", base=2)
ax1.set_yscale("log", base=10)
ax1.set_ylim([1e-16, 0.2e-15])
ax1.set_xlabel("$N$")
ax1.set_ylabel(r"$|ϵ|_2$")
ax1.set_title("Spectral method")

ax2.plot(grid, fft_fdm, color="blue",lw=4,
         marker = "o", markeredgecolor="k",
         markersize=12)
ax2.plot(ref_x, ref_y, color="black",lw=3, ls ="--")
ax2.set_xscale("log", base=2)
ax2.set_yscale("log", base=10)
ax2.set_ylim([4e-6, 4e-3])
ax2.set_xlabel("$N$")
ax2.set_ylabel(r"$|ϵ|_2$")
ax2.set_title("Second-order CDS")

fig.tight_layout()
fig.savefig("order_FFT.pdf")
plt.show()

