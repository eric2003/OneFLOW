import numpy as np
import matplotlib.pyplot as plt

# 读取数据
data = np.loadtxt('precise_wave.dat', comments='#')
x = data[:, 0]
u = data[:, 1]

# 简单绘图
plt.figure(figsize=(10, 6))
plt.plot(x, u, 'b-', linewidth=2)
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Complex Wave')
plt.grid(True)
plt.show()