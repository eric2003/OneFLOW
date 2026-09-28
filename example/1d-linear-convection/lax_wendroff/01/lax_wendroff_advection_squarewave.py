import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置（与之前一致，便于对比）
nx = 101                 # 空间网格点数
L = 4.0                  # 空间区间长度 [0, L]
x = np.linspace(0, L, nx)
dx = x[1] - x[0]         # dx = 0.04
c = 1.0                  # 对流速度（c > 0，向右传播）
dt = 0.01                # 时间步长（CFL ≈ 0.25，满足 Lax-Wendroff 稳定性 CFL ≤ 1）
total_time = 1.0         # 总模拟时间（位移 = 1.0）
nt = int(total_time / dt + 1e-8)   # 时间步数 ≈ 100

# 初始条件：方波
def initial_square_wave(x):
    u = np.zeros_like(x)
    u[(x >= 0.5) & (x <= 1.5)] = 1.0
    return u

u_initial = initial_square_wave(x)

# Lax-Wendroff 中心格式推进函数（显式、二阶时间二阶空间、矢量化实现）
def lax_wendroff(u_initial, dx, dt, c, nt):
    u = u_initial.copy()
    courant = c * dt / dx
    alpha = courant**2 / 2
    beta = courant / 2
    for n in range(nt):
        u_ip1 = np.roll(u, -1)
        u_im1 = np.roll(u, 1)
        u = u - beta * (u_ip1 - u_im1) + alpha * (u_ip1 - 2*u + u_im1)
    return u

# 计算解
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))   # 理论解（精确平移）
u_central = lax_wendroff(u_initial, dx, dt, c, nt)             # Lax-Wendroff 中心格式

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

print_errors("Lax-Wendroff 中心格式数值解", u_central, u_exact, dx)

# 绘图：理论解黑色实线，中心格式绿色三角符号（带细连线）
plt.figure(figsize=(12, 7))
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解 (精确平移)')
plt.plot(x, u_central, 'g^-', lw=1, markersize=6, markevery=2, label='Lax-Wendroff 中心格式数值解')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程 Lax-Wendroff 中心格式：方波传播最终对比 (nx=101, t=1.0)')
plt.legend(fontsize=12)
plt.grid(True)
plt.ylim(-0.3, 1.4)   # 稍放大以显示振荡
plt.show()