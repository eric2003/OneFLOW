import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置（与之前一致）
nx = 401                # 空间网格点数
x = np.linspace(0, 4.0, nx)   # 空间区间 [0, 4]
dx = x[1] - x[0]         # dx ≈ 0.01
c = 1.0                  # 对流速度（c > 0，向右传播）
dt = 0.005               # 时间步长（CFL ≈ 0.5）
total_time = 1.0         # 总模拟时间
nt = int(total_time / dt)   # 时间步数 = 200

# 初始条件：方波
u = np.zeros(nx)
u[(x >= 0.5) & (x <= 1.5)] = 1.0
u_initial = u.copy()

# minmod 限幅器（矢量化）
def minmod(a, b):
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))

# 时间推进：二阶迎风 TVD 格式（minmod 限幅 + 显式 Euler）
for n in range(nt):
    # 计算左右斜率
    delta_left = (u - np.roll(u, 1)) / dx      # (u_i - u_{i-1}) / dx
    delta_right = (np.roll(u, -1) - u) / dx    # (u_{i+1} - u_i) / dx
    sigma = minmod(delta_left, delta_right)   # 限幅斜率
    
    # 左侧外推到界面 i+1/2 的值（c > 0，使用左迎风）
    u_interface = u + 0.5 * sigma * dx
    
    # 通量（flux = c * u）
    flux = c * u_interface
    
    # 更新（矢量化）
    u = u - (c * dt / dx) * (flux - np.roll(flux, 1))

# 理论解：精确平移
displacement = c * total_time
shift = int(round(displacement / dx))
u_exact = np.roll(u_initial, shift)

# 计算并输出误差
error_linf = np.max(np.abs(u - u_exact))
error_l1   = np.sum(np.abs(u - u_exact)) * dx
error_l2   = np.sqrt(np.sum((u - u_exact)**2) * dx)

print("=== 数值解与理论解误差统计 (t=1.0) ===")
print(f"L∞ 误差 (最大绝对误差): {error_linf:.6f}")
print(f"L1 误差: {error_l1:.6f}")
print(f"L2 误差: {error_l2:.6f}")
print(f"数值解最大值: {np.max(u):.6f}, 最小值: {np.min(u):.6f}")
print(f"理论解最大值: {np.max(u_exact):.6f}, 最小值: {np.min(u_exact):.6f}")

# 绘图：理论解用黑色实线，数值解用红色圆圈符号（带连线，便于观察网格点值）
plt.figure(figsize=(10, 6))
plt.plot(x, u_exact, 'k-', lw=2, label='理论解 (精确平移)')          # 黑色实线
plt.plot(x, u, 'ro-', lw=1, markersize=4, markevery=5, label='二阶迎风（minmod 限幅）数值解 (t=1.0)')  # 红色圆圈 + 细连线，markevery=5 稀疏显示符号避免拥挤
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程 二阶迎风格式（带 minmod 限幅器）：方波传播最终对比 (t=1.0)')
plt.legend()
plt.grid(True)
plt.ylim(-0.2, 1.2)
plt.show()