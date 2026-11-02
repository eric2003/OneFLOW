import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# 创建图形和轴
fig, ax = plt.subplots()
ax.set_xlim(0, 2*np.pi)
ax.set_ylim(-1, 1)

# 初始化数据和时间标签
x = np.linspace(0, 2*np.pi, 1000)
line, = ax.plot(x, np.sin(x), color='blue')
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=12)

# 动画更新函数：引入时间变化
def animate(frame):
    t = frame * 0.1  # 时间以秒为单位（每帧 0.1 秒）
    y = np.sin(x + t)  # 波形随时间偏移
    line.set_ydata(y)
    time_text.set_text(f'时间: {t:.1f} s')  # 显示当前时间
    return line, time_text

# 创建动画：200 帧，每帧间隔 50ms（总时长约 10 秒）
ani = animation.FuncAnimation(fig, animate, frames=200, interval=50, blit=True)

# 显示动画
plt.show()

# 保存为 GIF（需要 Pillow: pip install pillow）
ani.save('sine_wave.gif', writer='pillow', fps=20)

# 保存为 MP4（需要 FFmpeg: 安装后配置环境变量）
ani.save('sine_wave.mp4', writer='ffmpeg', fps=20)