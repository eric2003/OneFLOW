import matplotlib.pyplot as plt
import numpy as np

# 设置图形大小和样式
plt.figure(figsize=(12, 6))

# 绘制 (a) Left-side reconstruction
plt.subplot(1, 2, 1)  # 1行2列的第一个子图

# 绘制网格点 (红点)
x_points = np.array([-2, -1, 0, 1, 2, 3, 4, 5, 6])  # 对应 i = -2 到 N+3
y_points = np.zeros_like(x_points)  # 所有点在 y=0 线上

plt.plot(x_points, y_points, 'ro', markersize=10)  # 红点

# 标记特定点 (i=1 和其他点)
plt.text(-2, 0.1, '$i=-2$', fontsize=12, ha='center')
plt.text(-1, 0.1, '$i=-1$', fontsize=12, ha='center')
plt.text(0, 0.1, '$i=0$', fontsize=12, ha='center')
plt.text(1, 0.1, '$i=1$', fontsize=12, ha='center')
plt.text(2, 0.1, '$i=2$', fontsize=12, ha='center')
plt.text(3, 0.1, '$i=3$', fontsize=12, ha='center')
plt.text(4, 0.1, '$i=N+1$', fontsize=12, ha='center')
plt.text(5, 0.1, '$i=N+2$', fontsize=12, ha='center')
plt.text(6, 0.1, '$i=N+3$', fontsize=12, ha='center')

# 绘制 x=0 和 x=L 的虚线
plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
plt.axvline(x=6, color='k', linestyle='--', alpha=0.5)
plt.text(0, -0.2, '$x=0$', fontsize=12, ha='center')
plt.text(6, -0.2, '$x=L$', fontsize=12, ha='center')

# 绘制紫色矩形模板 (i-2 到 i+2, 覆盖 i+1)
template_left = plt.Rectangle((0.5, -0.5), 4, 1, color='purple', alpha=0.3)
plt.gca().add_patch(template_left)
plt.text(2.5, 0.3, '$u_{i+1}^L$', fontsize=12, ha='center')

# 设置标题和轴
plt.title('(a) Left-side reconstruction', fontsize=14)
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 绘制 (b) Right-side reconstruction
plt.subplot(1, 2, 2)  # 1行2列的第二个子图

# 绘制网格点 (红点)
plt.plot(x_points, y_points, 'ro', markersize=10)  # 红点

# 标记特定点 (i=1 和其他点)
plt.text(-2, 0.1, '$i=-2$', fontsize=12, ha='center')
plt.text(-1, 0.1, '$i=-1$', fontsize=12, ha='center')
plt.text(0, 0.1, '$i=0$', fontsize=12, ha='center')
plt.text(1, 0.1, '$i=1$', fontsize=12, ha='center')
plt.text(2, 0.1, '$i=2$', fontsize=12, ha='center')
plt.text(3, 0.1, '$i=3$', fontsize=12, ha='center')
plt.text(4, 0.1, '$i=N+1$', fontsize=12, ha='center')
plt.text(5, 0.1, '$i=N+2$', fontsize=12, ha='center')
plt.text(6, 0.1, '$i=N+3$', fontsize=12, ha='center')

# 绘制 x=0 和 x=L 的虚线
plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
plt.axvline(x=6, color='k', linestyle='--', alpha=0.5)
plt.text(0, -0.2, '$x=0$', fontsize=12, ha='center')
plt.text(6, -0.2, '$x=L$', fontsize=12, ha='center')

# 绘制紫色矩形模板 (i-1 到 i+3, 覆盖 i+1)
template_right = plt.Rectangle((0.5, -0.5), 5, 1, color='purple', alpha=0.3)
plt.gca().add_patch(template_right)
plt.text(2.5, 0.3, '$u_{i+1}^R$', fontsize=12, ha='center')

# 绘制黑点 (i=1)
plt.plot(1, 0, 'ko', markersize=10)  # 黑点在 i=1

# 设置标题和轴
plt.title('(b) Right-side reconstruction', fontsize=14)
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()