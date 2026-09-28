import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(8, 4))  # 长方形画布

# 画十字线
plt.plot([-6, 6], [0, 0], 'r-', linewidth=2)  # 水平线
plt.plot([0, 0], [-2, 2], 'b-', linewidth=2)    # 垂直线

# 设置坐标范围
plt.xlim(-6, 6)
plt.ylim(-2, 2)

# 关键：使用默认的宽高比，不强制等比例
ax = plt.gca()
ax.set_aspect('auto')  # 这是默认设置，可以省略

plt.title('Auto aspect ratio (default)')
plt.grid(True)
plt.show()