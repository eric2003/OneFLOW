import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

# 1. 创建主窗口
root = tk.Tk()
root.title("Matplotlib嵌入Tkinter示例")
root.geometry("800x600")  # 设置窗口大小

# 2. 创建Matplotlib的绘图对象
# 创建一个Figure实例，指定尺寸和分辨率
fig = Figure(figsize=(7, 5), dpi=100)
# 添加子图（1行1列第1个）
ax = fig.add_subplot(111)

# 3. 绘制示例图形（正弦曲线）
x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)
ax.plot(x, y, label='sin(x)', color='blue', linewidth=2)
ax.set_title("正弦曲线", fontsize=14)
ax.set_xlabel("X轴", fontsize=12)
ax.set_ylabel("Y轴", fontsize=12)
ax.legend()
ax.grid(True, alpha=0.3)

# 4. 将Matplotlib图形嵌入Tkinter窗口
# 创建画布对象，关联Figure和Tkinter窗口
canvas = FigureCanvasTkAgg(fig, master=root)
# 绘制画布
canvas.draw()
# 将画布包装成Tkinter的组件并布局
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# 5. 可选：添加界面交互按钮（示例：重新绘制余弦曲线）
def redraw_cos():
    ax.clear()  # 清空原有图形
    x = np.linspace(0, 2 * np.pi, 100)
    y = np.cos(x)
    ax.plot(x, y, label='cos(x)', color='red', linewidth=2)
    ax.set_title("余弦曲线", fontsize=14)
    ax.set_xlabel("X轴", fontsize=12)
    ax.set_ylabel("Y轴", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    canvas.draw()  # 重新绘制画布

# 添加按钮
btn = tk.Button(root, text="切换为余弦曲线", command=redraw_cos)
btn.pack(side=tk.BOTTOM, pady=10)

# 6. 运行主循环
if __name__ == "__main__":
    root.mainloop()