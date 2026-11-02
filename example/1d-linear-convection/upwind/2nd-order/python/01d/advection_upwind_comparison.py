import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置（网格点减少到 101，便于快速运行和观察）
nx = 101                 # 空间网格点数
L = 4.0                  # 空间区间长度 [0, L]
x = np.linspace(0, L, nx)
dx = x[1] - x[0]         # dx = 0.04
c = 1.0                  # 对流速度（c > 0，向右传播）
dt = 0.01                # 时间步长（CFL ≈ 0.25，保守稳定）
total_time = 1.0         # 总模拟时间（位移 = 1.0）
nt = int(total_time / dt + 1e-8)   # 时间步数 ≈ 100

# 初始条件：方波
def initial_square_wave(x):
    u = np.zeros_like(x)
    u[(x >= 0.5) & (x <= 1.5)] = 1.0
    return u

u_initial = initial_square_wave(x)

# minmod 限幅器（矢量化）
def minmod(a, b):
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))

# 一阶迎风推进函数（显式、矢量化实现）
def first_order_upwind(u_initial, dx, dt, c, nt):
    u = u_initial.copy()
    courant = c * dt / dx
    for n in range(nt):
        u = u - courant * (u - np.roll(u, 1))   # 后向差分
    return u

# 二阶迎风 TVD（minmod 限幅）推进函数（显式、矢量化实现）
def second_order_tvd_minmod(u_initial, dx, dt, c, nt):
    u = u_initial.copy()
    for n in range(nt):
        delta_left  = (u - np.roll(u, 1)) / dx
        delta_right = (np.roll(u, -1) - u) / dx
        sigma = minmod(delta_left, delta_right)          # 限幅斜率
        u_interface = u + 0.5 * dx * sigma               # 左侧外推到 i+1/2 界面
        flux = c * u_interface
        u = u - (c * dt / dx) * (flux - np.roll(flux, 1))
    return u

# 计算三种解
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))   # 理论解（精确平移）

u1 = first_order_upwind(u_initial, dx, dt, c, nt)              # 一阶迎风
u2 = second_order_tvd_minmod(u_initial, dx, dt, c, nt)         # 二阶迎风（minmod TVD）

# 计算并输出误差
def print_errors(label, u_num, u_exact, dx):
    error_linf = np.max(np.abs(u_num - u_exact))
    error_l1   = np.sum(np.abs(u_num - u_exact)) * dx
    error_l2   = np.sqrt(np.sum((u_num - u_exact)**2) * dx)
    print(f"=== {label} 误差统计 (t={total_time}) ===")
    print(f"L∞ 误差: {error_linf:.6f}")
    print(f"L1 误差: {error_l1:.6f}")
    print(f"L2 误差: {error_l2:.6f}")
    print(f"最大值: {np.max(u_num):.6f}, 最小值: {np.min(u_num):.6f}\n")

print_errors("一阶迎风数值解", u1, u_exact, dx)
print_errors("二阶迎风（minmod TVD）数值解", u2, u_exact, dx)

# 绘图：理论解黑色实线，一阶蓝色圆圈，二阶红色方块（带细连线，便于区分）
plt.figure(figsize=(12, 7))
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解 (精确平移)')
plt.plot(x, u1, 'bo-', lw=1, markersize=6, markevery=2, label='一阶迎风数值解')
plt.plot(x, u2, 'rs-', lw=1, markersize=6, markevery=2, label='二阶迎风（minmod TVD）数值解')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程：1阶 vs 2阶迎风格式 方波传播对比 (nx=101, t=1.0)')
plt.legend(fontsize=12)
plt.grid(True)
plt.ylim(-0.2, 1.3)
plt.show()