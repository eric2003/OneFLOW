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

plt.show()

# 进度回调函数：实时显示保存进度
def progress_callback(i, n):
    print(f'保存进度: {i}/{n} 帧 ({i/n*100:.1f}%)')

# 先显示动画（非阻塞窗口）
plt.show(block=False)  # 窗口弹出播放，但脚本继续运行
print('动画窗口已打开，正在播放...')

# 再保存文件（带进度和错误处理）
try:
    # 选择保存格式：取消注释需要的
    ani.save('sine_wave.gif', writer='pillow', fps=20, progress_callback=progress_callback)   # 保存 GIF
    # ani.save('sine_wave.mp4', writer='ffmpeg', fps=20, progress_callback=progress_callback)  # 保存 MP4（需 FFmpeg）
    print('保存完毕！')
except Exception as e:
    print(f'保存出错: {e}')
    print('请检查 writer（如 Pillow/FFmpeg）是否安装，或帧数是否过大。')

# 保持窗口打开（可选，防止立即关闭）
plt.pause(0.01)  # 短暂暂停，让窗口持续