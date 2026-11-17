import matplotlib.pyplot as plt
import numpy as np

# 创建图形
fig, ax = plt.subplots(figsize=(6, 6))

# 绘制向上的箭头
# plt.arrow(x起点, y起点, dx水平增量, dy垂直增量, ...)
ax.arrow(0.5, 0.2, 0, 0.6, 
         head_width=0.05, head_length=0.1, 
         fc='blue', ec='blue', linewidth=2)

# 设置坐标轴范围
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')  # 等比例显示
plt.title('向上的箭头 - plt.arrow()')
plt.show()