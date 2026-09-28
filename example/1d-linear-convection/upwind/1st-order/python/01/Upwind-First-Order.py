import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置（与之前 Crank-Nicolson 保持一致，便于对比）
nx = 401                # 空间网格点数
x = np.linspace(0, 4.0, nx)   # 空间区间 [0, 4]
dx = x[1] - x[0]         # dx ≈ 0.01
c = 1.0                  # 对流速度（c > 0，向右传播，使用左迎风）
dt = 0.005               # 时间步长（CFL = c*dt/dx ≈ 0.5 < 1，满足稳定性）
total_time = 1.0         # 总模拟时间（位移 = 1.0）
nt = int(total_time / dt)   # 时间步数 = 200

courant = c * dt / dx     # Courant 数 ν ≈ 0.5

# 初始条件：方波
u = np.zeros(nx)
u[(x >= 0.5) & (x <= 1.5)] = 1.0
u_initial = u.copy()

# 时间推进：显式一阶迎风格式（c > 0，使用后向差分）
for n in range(nt):
    u_new = u.copy()
    for i in range(nx):
        im1 = (i - 1) % nx    # 周期边界
        u_new[i] = u[i] - courant * (u[i] - u[im1])
    u = u_new

# 理论解：精确平移
displacement = c * total_time
shift = int(round(displacement / dx))
u_exact = np.roll(u_initial, shift)

# 绘图：只显示最终时刻的数值解和理论解
plt.figure(figsize=(10, 6))
plt.plot(x, u, 'r-', lw=2, label='1阶迎风数值解 (t=1.0)')
plt.plot(x, u_exact, 'k--', lw=2, label='理论解 (精确平移)')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程 1阶迎风格式：方波传播最终对比 (t=1.0)')
plt.legend()
plt.grid(True)
plt.ylim(-0.2, 1.2)
plt.show()