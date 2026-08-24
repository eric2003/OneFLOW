import tkinter as tk
from tkinter import ttk, Checkbutton, StringVar, BooleanVar
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

class PlotGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("可交互绘图界面（线形/显隐控制）")
        self.root.geometry("900x700")

        # 1. 初始化数据和状态变量
        self.x = np.linspace(0, 2 * np.pi, 200)
        # 定义多条曲线的数据
        self.curves = {
            "正弦曲线": np.sin(self.x),
            "余弦曲线": np.cos(self.x),
            "正切曲线": np.tan(self.x)
        }
        # 曲线显隐状态（默认都显示）
        self.curve_visible = {name: BooleanVar(value=True) for name in self.curves.keys()}
        # 线形选择（默认实线）
        self.line_style = StringVar(value="-")
        # 存储绘制的曲线对象，方便后续修改样式/显隐
        self.line_objects = {}

        # 2. 创建界面布局
        self.create_widgets()
        # 3. 初始化绘图
        self.init_plot()
        # 4. 首次绘制曲线
        self.update_plot()

    def create_widgets(self):
        # ===== 左侧控制面板 =====
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        # 线形选择下拉框
        ttk.Label(control_frame, text="选择线形：").pack(anchor=tk.W, pady=5)
        line_styles = [("实线", "-"), ("虚线", "--"), ("点线", ":"), ("点划线", "-.")]
        style_combo = ttk.Combobox(
            control_frame, textvariable=self.line_style, values=[s[0] for s in line_styles], state="readonly"
        )
        style_combo.pack(anchor=tk.W, fill=tk.X, pady=5)
        style_combo.bind("<<ComboboxSelected>>", self.on_style_change)

        # 曲线显隐勾选框
        ttk.Label(control_frame, text="显示/隐藏曲线：").pack(anchor=tk.W, pady=10)
        for name in self.curves.keys():
            cb = Checkbutton(
                control_frame, text=name, variable=self.curve_visible[name],
                command=self.on_visibility_change
            )
            cb.pack(anchor=tk.W, pady=2)

        # ===== 右侧绘图区域 =====
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 创建Matplotlib绘图对象
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("X轴", fontsize=12)
        self.ax.set_ylabel("Y轴", fontsize=12)
        self.ax.set_title("可交互曲线控制", fontsize=14)
        self.ax.grid(True, alpha=0.3)
        self.ax.set_ylim(-3, 3)  # 限制y轴范围（避免正切曲线无限大）

        # 绑定到Tkinter画布
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def init_plot(self):
        """初始化曲线对象（先绘制所有曲线，后续控制显隐/样式）"""
        for name, y_data in self.curves.items():
            # 绘制曲线并保存对象引用
            line, = self.ax.plot(self.x, y_data, label=name, linestyle=self.line_style.get())
            self.line_objects[name] = line

    def update_plot(self):
        """更新所有曲线的显隐状态和样式"""
        for name, line in self.line_objects.items():
            # 设置显隐
            line.set_visible(self.curve_visible[name].get())
            # 设置线形
            line.set_linestyle(self.line_style.get())
        # 更新图例（同步显隐状态）
        self.ax.legend(loc="upper right")
        # 重新绘制画布
        self.canvas.draw()

    def on_style_change(self, event):
        """线形选择变化时的回调"""
        self.update_plot()

    def on_visibility_change(self):
        """勾选框变化时的回调（显隐曲线）"""
        self.update_plot()

if __name__ == "__main__":
    root = tk.Tk()
    app = PlotGUI(root)
    root.mainloop()