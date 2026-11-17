import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置
nx = 401                # 空间网格点数
x = np.linspace(0, 4.0, nx)   # 空间区间 [0, 4]
dx = x[1] - x[0]         # 空间步长 dx = 0.01
c = 1.0                  # 对流速度（正值向右传播）
dt = 0.005               # 时间步长（CFL ≈ 0.5）
total_time = 1.0         # 总模拟时间（位移 = c * t = 1.0，正好移动 100 个网格）
nt = int(total_time / dt)   # 时间步数 = 200

r = c * dt / (4 * dx)    # Crank-Nicolson 参数 r ≈ 0.125

# 初始条件：方波（x 在 [0.5, 1.5] 区间为 1，其余为 0）
u = np.zeros(nx)
u[(x >= 0.5) & (x <= 1.5)] = 1.0
u_initial = u.copy()

# 构建左端矩阵 A（固定不变，周期边界）
A = np.eye(nx)           # 对角线为 1
for i in range(nx):
    ip1 = (i + 1) % nx    # 周期边界
    im1 = (i - 1) % nx
    A[i, ip1] = r         # +r 项
    A[i, im1] = -r        # -r 项

# 时间循环
for n in range(nt):
    b = np.zeros(nx)
    for i in range(nx):
        ip1 = (i + 1) % nx
        im1 = (i - 1) % nx
        # 右端向量 b： r * u_{i-1} + u_i - r * u_{i+1}
        b[i] = r * u[im1] + u[i] - r * u[ip1]
    # 求解线性方程组 A u^{n+1} = b
    u = np.linalg.solve(A, b)

# 理论解：初始方波整体平移（位移 1.0 对应正好 100 个网格点）
displacement = c * total_time
shift = int(round(displacement / dx))   # 精确整数网格移动
u_exact = np.roll(u_initial, shift)

# 绘图对比
plt.figure(figsize=(10, 6))
plt.plot(x, u_initial, 'b-', lw=2, label='初始条件')
plt.plot(x, u, 'r-', lw=2, label='Crank-Nicolson 数值解 (t=1.0)')
plt.plot(x, u_exact, 'k--', lw=2, label='理论解 (精确平移)')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程 Crank-Nicolson 格式：方波传播对比 (t = %.1f)' % total_time)
plt.legend()
plt.grid(True)
plt.ylim(-0.2, 1.2)   # 放大 y 轴以显示可能的轻微振荡
plt.show()