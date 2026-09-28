import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# 创建图形和轴
fig, ax = plt.subplots()
ax.set_xlim(0, 2*np.pi)
ax.set_ylim(-1, 1)

# 初始化数据
x = np.linspace(0, 2*np.pi, 1000)
line, = ax.plot(x, np.sin(x), color='blue')

# 动画更新函数
def animate(frame):
    # 随着帧数增加，x 数据偏移，实现波形移动
    y = np.sin(x + frame / 10.0)
    line.set_ydata(y)
    return line,

# 创建动画：100 帧，每帧间隔 20ms，支持 blit 优化（更快渲染）
ani = animation.FuncAnimation(fig, animate, frames=100, interval=20, blit=True)

# 显示动画（在 Jupyter 中可用 plt.show()，否则保存为 GIF）
plt.show()

# 可选：保存为 GIF 文件（需要 pillow 或 imagemagick）
# ani.save('sine_wave.gif', writer='pillow', fps=30)