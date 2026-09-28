import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 参数：加密空间网格
nx = 1601                # 加密到 1601 点（dx ≈ 0.0025）
x = np.linspace(0, 4.0, nx)
dx = x[1] - x[0]
c = 1.0
dt = 0.005               # dt 不变，CFL ≈ 2（仍稳定）
total_time = 1.0
nt = int(total_time / dt)

r = c * dt / (4 * dx)

# 初始条件
u = np.zeros(nx)
u[(x >= 0.5) & (x <= 1.5)] = 1.0
u_initial = u.copy()

# 矩阵 A
A = np.eye(nx)
for i in range(nx):
    ip1 = (i + 1) % nx
    im1 = (i - 1) % nx
    A[i, ip1] = r
    A[i, im1] = -r

# 时间推进（注意：nx大后计算稍慢，但仍可接受）
for n in range(nt):
    b = np.zeros(nx)
    for i in range(nx):
        ip1 = (i + 1) % nx
        im1 = (i - 1) % nx
        b[i] = r * u[im1] + u[i] - r * u[ip1]
    u = np.linalg.solve(A, b)

# 理论解
shift = int(round(c * total_time / dx))
u_exact = np.roll(u_initial, shift)

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(x, u, 'r-', lw=2, label='Crank-Nicolson 数值解 (t=1.0)')
plt.plot(x, u_exact, 'k--', lw=2, label='理论解 (精确平移)')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程 Crank-Nicolson（加密网格 nx=1601）最终对比 (t=1.0)')
plt.legend()
plt.grid(True)
plt.ylim(-0.2, 1.2)
plt.show()