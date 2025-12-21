import matplotlib.pyplot as plt
import numpy as np

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# 左图：简单边界
# 主竖线
boundary_x = -1
ax1.axvline(x=boundary_x, color='black', linewidth=4)

# 添加斜线装饰
num_slashes = 12
for y in np.linspace(-2, 2, num_slashes):
    # 创建短斜线簇（每组3条）
    for offset in [-0.1, 0, 0.1]:
        ax1.plot([boundary_x, boundary_x + 0.5],
                 [y + offset, y + offset + 0.3],
                 color='darkred', linewidth=1.5, alpha=0.7)

ax1.set_xlim(-2, 1)
ax1.set_ylim(-3, 3)
ax1.set_aspect('equal')
ax1.set_title('Simple boundary with slash clusters')
ax1.grid(True, alpha=0.3)

# 右图：复杂边界模式
# 多条竖线形成边界带
for x in np.arange(-1, 1.1, 0.2):
    ax2.axvline(x=x, color='gray', linewidth=1, alpha=0.3)

# 主边界线
main_boundary = 0
ax2.axvline(x=main_boundary, color='darkblue', linewidth=3)

# 在边界上添加有规律的斜线
y_positions = np.arange(-2.5, 3, 0.5)
for i, y in enumerate(y_positions):
    # 交替的斜线方向
    direction = -1 if i % 2 == 0 else 1
    angle = 45 * direction
    
    # 画斜线
    length = 0.4
    x_end = main_boundary + length * np.cos(np.radians(angle))
    y_end = y + length * np.sin(np.radians(angle))
    
    ax2.plot([main_boundary, x_end], [y, y_end],
             color='orange', linewidth=2, alpha=0.8)
    
    # 在斜线末端添加箭头
    ax2.scatter(x_end, y_end, color='red', s=30, alpha=0.7)

# 添加文字说明
ax2.text(main_boundary + 0.5, 2.5, 'Field/Area', 
         fontsize=12, ha='center', va='center',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
ax2.text(main_boundary - 0.5, 2.5, 'Boundary', 
         fontsize=12, ha='center', va='center',
         bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))

ax2.set_xlim(-1.5, 1.5)
ax2.set_ylim(-3, 3)
ax2.set_aspect('equal')
ax2.set_title('Complex boundary pattern with directional indicators')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()