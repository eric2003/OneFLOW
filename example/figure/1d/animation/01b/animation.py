import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# 设置中文字体（Windows 推荐）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']  # 支持中文
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 创建图形和轴
fig, ax = plt.subplots()
ax.set_xlim(0, 2*np.pi)
ax.set_ylim(-1, 1)
ax.set_title('正弦波动画')  # 测试中文标题

# 初始化数据和时间标签
x = np.linspace(0, 2*np.pi, 1000)
line, = ax.plot(x, np.sin(x), color='blue')
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=12)

# 动画更新函数：引入时间变化
def animate(frame):
    t = frame * 0.1  # 时间以秒为单位
    y = np.sin(x + t)
    line.set_ydata(y)
    time_text.set_text(f'时间: {t:.1f} s')  # 中文现在正常
    return line, time_text

# 创建动画
ani = animation.FuncAnimation(fig, animate, frames=200, interval=50, blit=True)

# 先保存文件（取消注释需要的）
ani.save('sine_wave.gif', writer='pillow', fps=20)   # 保存 GIF
ani.save('sine_wave.mp4', writer='ffmpeg', fps=20)  # 保存 MP4（需 FFmpeg）

# 再显示动画（窗口播放）
plt.show()