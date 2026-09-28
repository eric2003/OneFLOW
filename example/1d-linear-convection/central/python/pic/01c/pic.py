import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 1. 生成数据
x = np.linspace(0, 2*np.pi, 200)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.tan(x)

# 2. 创建绘图窗口，开启交互模式
plt.ion()  # 开启交互模式
fig, ax = plt.subplots(figsize=(8, 6))

# 3. 绘制曲线并保存对象（用于后续控制）
line1, = ax.plot(x, y1, label='正弦曲线', picker=True)  # picker=True开启选中功能
line2, = ax.plot(x, y2, label='余弦曲线', picker=True)
line3, = ax.plot(x, y3, label='正切曲线', picker=True)

# 4. 基础配置
ax.set_title("Matplotlib原生交互示例", fontsize=14)
ax.set_ylim(-3, 3)
ax.legend()
ax.grid(alpha=0.3)

# ===== 实现核心交互功能 =====
def on_click(event):
    """点击图例/曲线，切换显隐"""
    if event.inaxes is None:
        return
    # 点击曲线切换显隐
    for line in [line1, line2, line3]:
        if line.contains(event)[0]:
            line.set_visible(not line.get_visible())
            plt.draw()
            return

def on_key(event):
    """键盘快捷键实现功能：改标题/导出数据/导出图片"""
    # 按t：修改标题
    if event.key == 't':
        ax.set_title("修改后的标题：自定义曲线交互", fontsize=14)
    # 按s：保存数据（当前显示的曲线数据）
    elif event.key == 's':
        with open('curve_data.txt', 'w') as f:
            f.write("x,sin(x),cos(x),tan(x)\n")
            for i in range(len(x)):
                f.write(f"{x[i]},{y1[i]},{y2[i]},{y3[i]}\n")
        print("数据已保存到curve_data.txt")
    # 按p：保存图片
    elif event.key == 'p':
        fig.savefig('plot.png', dpi=150, bbox_inches='tight')
        print("图片已保存到plot.png")
    plt.draw()

# 绑定事件
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_key)

# 5. 保持窗口显示（缩放/平移用Matplotlib原生工具栏）
plt.ioff()  # 关闭交互模式，防止窗口闪退
plt.show()