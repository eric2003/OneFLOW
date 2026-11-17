import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置（nx 不变）
nx = 101
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]         # dx ≈ 0.04
c = 1.0
total_time = 1.0

# 改进1: CFL=1
dt = dx / c              # dt ≈ 0.04, CFL=1（减小色散）
nt = int(total_time / dt + 1e-8)

# 改进2: 人工粘性系数（可调，0.01~0.05，越大耗散越强但波形越模糊）
epsilon = 0.02           # 四阶人工粘性强度

# 初始条件
def initial_square_wave(x):
    u = np.zeros_like(x)
    u[(x >= 0.5) & (x <= 1.5)] = 1.0
    return u

u_initial = initial_square_wave(x)

# 改进版 Lax-Wendroff（添加四阶人工粘性）
def lax_wendroff_with_av(u_initial, dx, dt, c, nt, epsilon):
    u = u_initial.copy()
    courant = c * dt / dx    # =1
    alpha = courant**2 / 2   # =0.5
    beta = courant / 2       # =0.5
    gamma = epsilon / dx     # 粘性系数缩放
    for n in range(nt):
        u_ip1 = np.roll(u, -1)
        u_im1 = np.roll(u, 1)
        u_ip2 = np.roll(u, -2)
        u_im2 = np.roll(u, 2)
        # 原 Lax-Wendroff
        advective = -beta * (u_ip1 - u_im1) + alpha * (u_ip1 - 2*u + u_im1)
        # 四阶人工粘性（抑制振荡）
        dissipative = gamma * (u_ip2 - 4*u_ip1 + 6*u - 4*u_im1 + u_im2)
        u = u + advective + dissipative
    return u

# 计算
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))
u_central = lax_wendroff_with_av(u_initial, dx, dt, c, nt, epsilon)

# 误差输出
def print_errors(label, u_num, u_exact, dx):
    error_linf = np.max(np.abs(u_num - u_exact))
    error_l1   = np.sum(np.abs(u_num - u_exact)) * dx
    error_l2   = np.sqrt(np.sum((u_num - u_exact)**2) * dx)
    print(f"=== {label} 误差统计 (t={total_time}) ===")
    print(f"L∞ 误差: {error_linf:.6f}")
    print(f"L1 误差: {error_l1:.6f}")
    print(f"L2 误差: {error_l2:.6f}")
    print(f"最大值: {np.max(u_num):.6f}, 最小值: {np.min(u_num):.6f}\n")

print_errors("改进 Lax-Wendroff (CFL=1 + 人工粘性)", u_central, u_exact, dx)

# 绘图
plt.figure(figsize=(12, 7))
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解 (精确平移)')
plt.plot(x, u_central, 'g^-', lw=1, markersize=6, markevery=2, label='改进 Lax-Wendroff 数值解')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程 改进 Lax-Wendroff：方波传播最终对比 (nx=101, t=1.0)')
plt.legend(fontsize=12)
plt.grid(True)
plt.ylim(-0.3, 1.4)
plt.show()