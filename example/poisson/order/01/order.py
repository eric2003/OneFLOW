import numpy as np
import matplotlib.pyplot as plt

grid = [32,64,128,256,512]
fft_spectral = [1.339154672220572e-16, 1.342489800590083e-16, 1.3269947532677886e-16, \
                1.4514549677756118e-16 , 1.485803753784319e-16]
fft_fdm = [0.0015607100315532957 , 0.0005987381110678801 , 0.00014313734718665358, \
           3.549617203207291e-5, 8.865373334924762e-6]
ref_x = [64,256]
ref_y = [2e-3, 1.25e-4]

#plt.figure("OneFLOW-CFD Solver+Order", figsize=(10, 4), dpi=100)
fig, ax = plt.subplots(figsize=(10, 4), dpi=100)
fig.canvas.manager.set_window_title("OneFLOW-CFD Solver+Order")
fig.suptitle('test title', fontsize=12)
ax.set_title('Fruit supply by kind and color')
ax.plot( grid, fft_spectral, color="red", marker = "o" )
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=10)
ax.set_ylim([1e-16, 0.2e-15])
ax.set_title("Spectral method")

plt.tight_layout()
plt.show()
