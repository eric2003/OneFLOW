import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 左图：长条形网格
ax1 = axes[0]
# 画主要网格线
for x in np.arange(-6, 7, 1):
    ax1.axvline(x=x, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
for y in np.arange(-2, 3, 0.5):
    ax1.axhline(y=y, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)

# 画粗的坐标轴
ax1.axhline(y=0, color='black', linewidth=2)
ax1.axvline(x=0, color='black', linewidth=2)

ax1.set_xlim(-6, 6)
ax1.set_ylim(-2, 2)
ax1.set_title('1D-like grid (long strip)')
ax1.set_xlabel('x-axis')
ax1.set_ylabel('y-axis')
ax1.set_aspect('auto')  # 不强制等比例

# 右图：比较等比例的情况
ax2 = axes[1]
# 同样画网格线
for x in np.arange(-6, 7, 1):
    ax2.axvline(x=x, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
for y in np.arange(-2, 3, 0.5):
    ax2.axhline(y=y, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)

ax2.axhline(y=0, color='black', linewidth=2)
ax2.axvline(x=0, color='black', linewidth=2)

ax2.set_xlim(-6, 6)
ax2.set_ylim(-2, 2)
ax2.set_title('With axis("equal") - becomes square')
ax2.set_xlabel('x-axis')
ax2.set_ylabel('y-axis')
ax2.set_aspect('equal')  # 强制等比例

plt.tight_layout()
plt.show()