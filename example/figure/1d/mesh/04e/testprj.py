import matplotlib.pyplot as plt
import numpy as np

# 最简单的边界面实现
plt.figure(figsize=(8, 6))

# 1. 画主竖线
plt.axvline(x=0, color='black', linewidth=4)

# 2. 在旁边添加斜线
y_positions = np.arange(-4, 4.5, 0.5)  # 从-4到4，步长0.5

for y in y_positions:
    # 每个位置画一条斜线
    plt.plot([0, 0.5],  # x坐标：从0到0.5
             [y, y + 0.3],  # y坐标：向上倾斜
             color='red', linewidth=2)

# 3. 设置显示范围
plt.xlim(-2, 2)
plt.ylim(-5, 5)

# 4. 添加标题和网格
plt.title('Simple boundary: vertical line with diagonal slashes')
plt.grid(True, alpha=0.3)

plt.show()