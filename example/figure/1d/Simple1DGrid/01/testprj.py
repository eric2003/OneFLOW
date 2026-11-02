import numpy as np
import matplotlib.pyplot as plt  # 可选，用于可视化

class Simple1DGrid:
    def __init__(self, a, b, N_physical, ng=2):
        """
        初始化1D网格。
        - a, b: 物理域边界 [x=a, x=b]
        - N_physical: 物理网格点数 (不含ghost)
        - ng: ghost层数 (每侧)
        """
        self.a = a
        self.b = b
        self.N = N_physical  # 物理点数
        self.ng = ng
        self.total_N = N_physical + 2 * ng  # 总数组大小
        
        # 网格点坐标 x (node-centered, 包括ghost)
        dx = (b - a) / N_physical
        self.x = np.linspace(a - ng * dx, b + ng * dx, self.total_N)
        
        # 网格中心坐标 xcc (cell-centered, 用于存储u等)
        self.xcc = 0.5 * (self.x[:-1] + self.x[1:])  # 总长 total_N - 1
        
        # 物理域切片 (视图，避免拷贝)
        self.ist = ng  # 逻辑起始 = ng (数组中物理从ng开始)
        self.ied = ng + N_physical  # 逻辑结束 = ng + N
        self.x_physical = self.x[self.ist:self.ied]
        self.xcc_physical = self.xcc[self.ist-1:self.ied-1]  # xcc物理部分 (N-1? 视方法而定)
        
        # 示例物理量: u (速度场，初始化为sin波)
        self.u = np.sin(2 * np.pi * self.xcc / (b - a))  # 全域初始化 (包括ghost，稍后填充)
        
        # 只暴露物理视图
        self.u_physical = self.u[self.ist-1:self.ied-1]  # 假设cell-centered，调整为self.ist:self.ied if node-centered
    
    def fill_ghosts_dirichlet(self, u_left=0.0, u_right=0.0):
        """
        填充ghost: Dirichlet边界 (u(a)=u_left, u(b)=u_right)
        """
        # 左侧ghost: 镜像或常值 (这里用常值)
        self.u[:self.ng] = u_left  # 简单常值填充
        # 右侧ghost
        self.u[-self.ng:] = u_right
        # 更新物理视图 (可选，但视图自动同步)
    
    def fill_ghosts_neumann(self, du_left=0.0, du_right=0.0):
        """
        Neumann边界 (du/dx = const)
        """
        dx = self.x[1] - self.x[0]
        # 左侧: u[-1] = u[0] - du_left * dx (一阶后向)
        self.u[:self.ng] = self.u[self.ng] - du_left * dx * np.arange(1, self.ng+1)[::-1]
        # 右侧类似
        self.u[-self.ng:] = self.u[-self.ng-1] + du_right * dx * np.arange(1, self.ng+1)
    
    def plot_grid(self):
        """可视化网格和u"""
        fig, ax = plt.subplots()
        ax.plot(self.x, np.zeros_like(self.x), 'ko-', label='All points (with ghosts)')
        ax.plot(self.x_physical, np.zeros_like(self.x_physical), 'ro-', label='Physical domain')
        ax.plot(self.xcc_physical, self.u_physical, 'b-', label='u at centers')
        ax.axvline(self.a, color='g', ls='--', label='x=a')
        ax.axvline(self.b, color='g', ls='--', label='x=b')
        ax.set_xlabel('x')
        ax.set_title('1D Grid with Ghosts')
        ax.legend()
        plt.show()

# 使用示例
if __name__ == "__main__":
    grid = Simple1DGrid(a=0.0, b=1.0, N_physical=50, ng=2)
    print(f"ist (offset): {grid.ist}, ied: {grid.ied}")
    print(f"Total array size: {grid.total_N}")
    print(f"Physical x range: {grid.x_physical[0]:.3f} to {grid.x_physical[-1]:.3f}")
    
    # 填充ghost (Dirichlet u=0)
    grid.fill_ghosts_dirichlet(u_left=0.0, u_right=0.0)
    
    # CFD时间步示例: 简单Advection (du/dt + c du/dx =0)，用upwind
    c = 0.5  # 波速
    dt = 0.01  # 时间步
    u_new = grid.u.copy()
    for i in range(grid.ist, grid.ied-1):  # 只循环物理域 (注意ied-1 for cell-centered)
        u_new[i] = grid.u[i] - c * dt * (grid.u[i] - grid.u[i-1]) / (grid.x[i] - grid.x[i-1])
    grid.u = u_new  # 更新
    
    grid.plot_grid()  # 可视化