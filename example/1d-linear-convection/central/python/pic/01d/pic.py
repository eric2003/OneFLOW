import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class SimpleTecplot:
    def __init__(self, root):
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号        
        self.root = root
        self.root.title("简易Tecplot曲线工具")
        self.root.geometry("1000x700")

        # 核心状态变量
        self.data = None  # 存储导入的数据
        self.curves = {}  # 存储曲线对象：{列名: line对象}
        self.visible = {} # 曲线显隐状态：{列名: BooleanVar}

        # 1. 创建界面（分面板：控制区+绘图区）
        self.create_ui()

        # 2. 初始化Matplotlib绘图
        self.init_plot()

    def create_ui(self):
        # ===== 左侧控制区 =====
        control_frame = ttk.Frame(self.root)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)  # 用padx/pady实现外间距

        # 1.1 数据导入按钮
        import_btn = ttk.Button(control_frame, text="导入数据(TXT/CSV/Excel)", 
                                command=self.import_data)
        import_btn.pack(fill=tk.X, pady=5)

        # 1.2 曲线显隐勾选框（动态生成）
        self.curve_frame = ttk.LabelFrame(control_frame, text="曲线控制")
        self.curve_frame.pack(fill=tk.BOTH, expand=True, pady=10)
        # 内部Frame实现内边距
        self.curve_inner_frame = ttk.Frame(self.curve_frame)
        self.curve_inner_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 1.3 样式调整面板
        style_frame = ttk.LabelFrame(control_frame, text="样式调整")
        style_frame.pack(fill=tk.X, pady=5)
        # 内部Frame实现内边距
        style_inner_frame = ttk.Frame(style_frame)
        style_inner_frame.pack(fill=tk.X, padx=5, pady=5)

        # 线型选择
        ttk.Label(style_inner_frame, text="线型：").pack(anchor=tk.W)
        self.linestyle = ttk.Combobox(style_inner_frame, values=["-", "--", ":", "-."])
        self.linestyle.set("-")
        self.linestyle.pack(fill=tk.X, pady=2)
        self.linestyle.bind("<<ComboboxSelected>>", self.update_style)

        # 线宽选择（最终修复：事件绑定用<ButtonRelease-1>）
        ttk.Label(style_inner_frame, text="线宽：").pack(anchor=tk.W)
        self.linewidth = tk.Scale(style_inner_frame, from_=0.5, to=5, orient=tk.HORIZONTAL)
        self.linewidth.set(1)  # 新版Python设置初始值的正确方式
        self.linewidth.pack(fill=tk.X, pady=2)
        # 修复事件绑定：Scale用<ButtonRelease-1>替代<Release>
        self.linewidth.bind("<ButtonRelease-1>", self.update_style)

        # 1.4 导出面板
        export_frame = ttk.LabelFrame(control_frame, text="导出")
        export_frame.pack(fill=tk.X, pady=5)
        # 内部Frame实现内边距
        export_inner_frame = ttk.Frame(export_frame)
        export_inner_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(export_inner_frame, text="导出图片", command=self.export_plot).pack(fill=tk.X, pady=2)
        ttk.Button(export_inner_frame, text="导出数据", command=self.export_data).pack(fill=tk.X, pady=2)

        # ===== 右侧绘图区 =====
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 创建Matplotlib画布
        self.fig, self.ax = plt.subplots(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def init_plot(self):
        """初始化绘图面板"""
        self.ax.set_xlabel("X轴")
        self.ax.set_ylabel("Y轴")
        self.ax.set_title("简易Tecplot曲线工具")
        self.ax.grid(alpha=0.3)
        self.canvas.draw()

    def import_data(self):
        """导入TXT/CSV/Excel数据"""
        # 选择文件
        file_path = filedialog.askopenfilename(
            filetypes=[("数据文件", "*.txt *.csv *.xlsx *.xls"), ("所有文件", "*.*")]
        )
        if not file_path:
            return

        # 读取数据（自动适配格式）
        try:
            if file_path.endswith((".txt", ".csv")):
                self.data = pd.read_csv(file_path, sep=None, engine="python")  # 自动识别分隔符
            elif file_path.endswith((".xlsx", ".xls")):
                self.data = pd.read_excel(file_path)
        except Exception as e:
            messagebox.showerror("错误", f"读取数据失败：{str(e)}")
            return

        # 清空原有曲线
        self.curves.clear()
        self.visible.clear()
        for widget in self.curve_inner_frame.winfo_children():
            widget.destroy()

        # 绘制曲线（默认用第一列做X轴，其他列做Y轴）
        x_col = self.data.columns[0]
        x = self.data[x_col]

        for y_col in self.data.columns[1:]:
            y = self.data[y_col]
            # 绘制曲线并保存对象
            line, = self.ax.plot(x, y, label=y_col, linestyle="-", linewidth=1)
            self.curves[y_col] = line
            # 创建显隐勾选框
            self.visible[y_col] = tk.BooleanVar(value=True)
            cb = tk.Checkbutton(self.curve_inner_frame, text=y_col, variable=self.visible[y_col],
                                 command=self.update_visibility)
            cb.pack(anchor=tk.W, pady=1)

        # 更新图例和画布
        self.ax.legend()
        self.canvas.draw()

    def update_visibility(self):
        """更新曲线显隐"""
        for col, var in self.visible.items():
            self.curves[col].set_visible(var.get())
        self.canvas.draw()

    def update_style(self, event=None):
        """更新曲线样式（线型/线宽）- 兼容无event的调用"""
        ls = self.linestyle.get()
        lw = self.linewidth.get()  # 获取Scale当前值
        for line in self.curves.values():
            line.set_linestyle(ls)
            line.set_linewidth(lw)
        self.canvas.draw()

    def export_plot(self):
        """导出图片"""
        save_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG图片", "*.png"), ("JPG图片", "*.jpg"), ("PDF文件", "*.pdf")]
        )
        if save_path:
            self.fig.savefig(save_path, dpi=150, bbox_inches="tight")

    def export_data(self):
        """导出当前显示的曲线数据"""
        if self.data is None:
            messagebox.showwarning("提示", "请先导入数据！")
            return
        save_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV文件", "*.csv"), ("TXT文件", "*.txt")]
        )
        if save_path:
            # 只导出显隐为True的列
            visible_cols = [self.data.columns[0]] + [col for col, var in self.visible.items() if var.get()]
            self.data[visible_cols].to_csv(save_path, index=False)

if __name__ == "__main__":
    root = tk.Tk()
    app = SimpleTecplot(root)
    root.mainloop()