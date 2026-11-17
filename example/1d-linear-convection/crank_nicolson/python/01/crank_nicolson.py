import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import time

class CrankNicolson1D:
    """
    一维线性对流方程的Crank-Nicolson格式求解器
    方程: ∂u/∂t + a ∂u/∂x = 0
    """
    
    def __init__(self, a=1.0, L=1.0, nx=101, cfl=0.5):
        """
        初始化参数
        
        参数:
        a: 对流速度 (默认 1.0)
        L: 计算域长度 (默认 1.0)
        nx: 空间网格数 (默认 101)
        cfl: CFL数 (默认 0.5，用于确定时间步长)
        """
        self.a = a  # 对流速度
        self.L = L  # 域长度
        self.nx = nx  # 空间网格数
        self.dx = L / (nx - 1)  # 空间步长
        self.x = np.linspace(0, L, nx)  # 空间网格
        
        # 根据CFL条件计算时间步长
        self.dt = cfl * self.dx / abs(a) if a != 0 else 0.01
        self.cfl = abs(a) * self.dt / self.dx
        
        print(f"空间步长 dx = {self.dx:.4f}")
        print(f"时间步长 dt = {self.dt:.6f}")
        print(f"CFL数 = {self.cfl:.4f}")
        
        # 初始化解
        self.u = np.zeros(nx)
        
    def set_initial_condition(self, ic_type='step'):
        """
        设置初始条件
        
        参数:
        ic_type: 初始条件类型
            'step': 阶跃函数
            'sine': 正弦波
            'gaussian': 高斯波包
        """
        if ic_type == 'step':
            # 阶跃函数 (Riemann问题)
            # u(x,0) = 1 if x < 0.5 else 0
            self.u = np.where(self.x < 0.5, 1.0, 0.0)
            
        elif ic_type == 'sine':
            # 正弦波
            self.u = np.sin(2 * np.pi * self.x)
            
        elif ic_type == 'gaussian':
            # 高斯波包
            x0 = 0.3
            sigma = 0.1
            self.u = np.exp(-((self.x - x0) / sigma) ** 2)
            
        elif ic_type == 'double_step':
            # 双阶跃函数，测试性更好
            self.u = np.zeros_like(self.x)
            self.u[(self.x >= 0.2) & (self.x <= 0.4)] = 1.0
            self.u[(self.x >= 0.6) & (self.x <= 0.8)] = 0.5
            
        else:
            raise ValueError(f"未知的初始条件类型: {ic_type}")
        
        # 保存初始条件用于对比
        self.u0 = self.u.copy()
        
    def assemble_matrix(self, boundary='periodic'):
        """
        组装Crank-Nicolson方法的线性系统矩阵
        
        参数:
        boundary: 边界条件类型
            'periodic': 周期性边界
            'dirichlet': Dirichlet边界
        """
        n = self.nx
        alpha = self.a * self.dt / (4 * self.dx)  # CN系数: 0.5 * dt * a / (2dx)
        
        # 主对角线: 1
        main_diag = np.ones(n)
        
        # 上对角线: -alpha (中心差分的贡献)
        upper_diag = -alpha * np.ones(n-1)
        
        # 下对角线: alpha
        lower_diag = alpha * np.ones(n-1)
        
        if boundary == 'periodic':
            # 周期性边界：需要额外的两个元素
            # A(0, n-1) 和 A(n-1, 0)
            diags_pos = [0, 1, -1, -(n-1), n-1]
            diags_vals = [main_diag, upper_diag, lower_diag, 
                         np.array([-alpha]), np.array([alpha])]
            
        elif boundary == 'dirichlet':
            # Dirichlet边界：边界点固定
            # 实际上对于对流方程，需要指定入流边界
            if self.a > 0:  # 向右传播
                # 左边界为入流边界，右边界为出流边界
                diags_pos = [0, 1, -1]
                diags_vals = [main_diag, upper_diag, lower_diag]
            else:  # 向左传播
                diags_pos = [0, 1, -1]
                diags_vals = [main_diag, upper_diag, lower_diag]
        else:
            raise ValueError(f"未知的边界条件: {boundary}")
        
        # 创建稀疏矩阵
        A = diags(diags_vals, diags_pos, shape=(n, n), format='csr')
        
        # 处理边界条件
        if boundary == 'dirichlet':
            if self.a > 0:
                # 左边界固定（入流边界）
                A[0, :] = 0
                A[0, 0] = 1
                # 右边界为自然边界（不做特殊处理）
            else:
                # 右边界固定（入流边界）
                A[-1, :] = 0
                A[-1, -1] = 1
        
        self.A = A
        return A
    
    def compute_rhs(self, u_old, boundary='periodic'):
        """
        计算右端项: b = [I + 0.5*dt*L] * u_old
        
        参数:
        u_old: 旧时间步的解
        boundary: 边界条件类型
        """
        n = self.nx
        alpha = self.a * self.dt / (4 * self.dx)  # 与矩阵组装相同的系数
        
        # 显式部分的计算
        b = u_old.copy()
        
        if boundary == 'periodic':
            # 周期性边界
            for i in range(n):
                # 中心差分
                left = u_old[(i-1) % n]
                right = u_old[(i+1) % n]
                b[i] += alpha * (right - left)
                
        elif boundary == 'dirichlet':
            # Dirichlet边界
            if self.a > 0:
                # 左边界固定
                b[0] = u_old[0]  # 保持边界值不变
                # 内点
                for i in range(1, n-1):
                    b[i] += alpha * (u_old[i+1] - u_old[i-1])
                # 右边界（出流）使用单边差分
                b[-1] += alpha * (u_old[-1] - u_old[-2])
            else:
                # 右边界固定
                b[-1] = u_old[-1]  # 保持边界值不变
                # 内点
                for i in range(1, n-1):
                    b[i] += alpha * (u_old[i+1] - u_old[i-1])
                # 左边界（出流）使用单边差分
                b[0] += alpha * (u_old[1] - u_old[0])
        
        return b
    
    def solve_step(self, u_old, boundary='periodic'):
        """
        求解一个时间步: A * u_new = b
        
        参数:
        u_old: 旧时间步的解
        boundary: 边界条件类型
        
        返回:
        u_new: 新时间步的解
        """
        # 组装矩阵（如果是第一次调用）
        if not hasattr(self, 'A'):
            self.assemble_matrix(boundary)
        
        # 计算右端项
        b = self.compute_rhs(u_old, boundary)
        
        # 处理边界条件的右端项
        if boundary == 'dirichlet':
            if self.a > 0:
                # 左边界固定值（假设为初始边界值）
                b[0] = self.u0[0]
            else:
                # 右边界固定值
                b[-1] = self.u0[-1]
        
        # 求解线性系统
        u_new = spsolve(self.A, b)
        
        return u_new
    
    def solve(self, T, boundary='periodic', ic_type='step', verbose=True):
        """
        从t=0到t=T求解
        
        参数:
        T: 总时间
        boundary: 边界条件类型
        ic_type: 初始条件类型
        verbose: 是否显示进度
        
        返回:
        u_final: 最终时刻的解
        history: 解的演化历史（可选）
        """
        # 设置初始条件
        self.set_initial_condition(ic_type)
        u = self.u.copy()
        
        # 计算时间步数
        nt = int(np.ceil(T / self.dt))
        actual_dt = T / nt  # 调整时间步长以恰好达到T
        
        if verbose:
            print(f"总时间 T = {T}")
            print(f"时间步数 = {nt}")
            print(f"调整后的时间步长 = {actual_dt:.6f}")
            print(f"实际CFL数 = {abs(self.a) * actual_dt / self.dx:.4f}")
        
        # 保存历史用于动画（每10步保存一次）
        save_interval = max(1, nt // 20)
        history = []
        history.append((0.0, u.copy()))
        
        # 时间推进
        start_time = time.time()
        
        for n in range(nt):
            # 使用实际的时间步长
            dt_actual = actual_dt if n == nt-1 else self.dt
            
            # 临时调整时间步长
            original_dt = self.dt
            self.dt = dt_actual
            
            # 重新组装矩阵（因为dt改变了）
            self.assemble_matrix(boundary)
            
            # 求解一个时间步
            u_new = self.solve_step(u, boundary)
            
            # 恢复原始时间步长
            self.dt = original_dt
            
            # 更新解
            u = u_new.copy()
            
            # 保存历史
            if (n+1) % save_interval == 0 or n == nt-1:
                current_time = (n+1) * dt_actual
                history.append((current_time, u.copy()))
                
                if verbose and (n+1) % (5*save_interval) == 0:
                    print(f"时间步 {n+1}/{nt}, t = {current_time:.3f}")
        
        end_time = time.time()
        
        if verbose:
            print(f"计算完成! 用时: {end_time - start_time:.2f} 秒")
        
        # 保存最终解
        self.u = u.copy()
        self.history = history
        
        return u, history
    
    def exact_solution(self, t, ic_type='step'):
        """
        计算精确解（对线性对流方程）
        
        参数:
        t: 时间
        ic_type: 初始条件类型
        
        返回:
        exact: 精确解
        """
        # 线性对流方程的精确解: u(x,t) = u0(x - a*t)
        x_exact = (self.x - self.a * t) % self.L  # 周期性边界
        
        # 根据初始条件类型计算精确解
        if ic_type == 'step':
            # 阶跃函数的精确解（周期性）
            exact = np.zeros_like(self.x)
            # 计算偏移后的阶跃位置
            x_shifted = (self.x - self.a * t) % self.L
            exact[x_shifted < 0.5] = 1.0
            
        elif ic_type == 'sine':
            # 正弦波的精确解
            exact = np.sin(2 * np.pi * x_exact)
            
        elif ic_type == 'gaussian':
            # 高斯波包的精确解
            x0 = 0.3
            sigma = 0.1
            x_shifted = (self.x - self.a * t) % self.L
            exact = np.exp(-((x_shifted - x0) / sigma) ** 2)
            
        elif ic_type == 'double_step':
            # 双阶跃的精确解
            exact = np.zeros_like(self.x)
            x_shifted = (self.x - self.a * t) % self.L
            exact[(x_shifted >= 0.2) & (x_shifted <= 0.4)] = 1.0
            exact[(x_shifted >= 0.6) & (x_shifted <= 0.8)] = 0.5
            
        return exact
    
    def compute_error(self, t, u_numeric):
        """
        计算数值解与精确解的误差
        
        参数:
        t: 时间
        u_numeric: 数值解
        
        返回:
        error_L1: L1误差
        error_L2: L2误差
        error_Linf: 最大误差
        """
        u_exact = self.exact_solution(t, self.ic_type if hasattr(self, 'ic_type') else 'step')
        
        error = u_numeric - u_exact
        error_L1 = np.sum(np.abs(error)) * self.dx
        error_L2 = np.sqrt(np.sum(error**2) * self.dx)
        error_Linf = np.max(np.abs(error))
        
        return error_L1, error_L2, error_Linf


def plot_results(solver, T=0.5, save_fig=False):
    """
    绘制结果并对比数值解和精确解
    """
    # 求解
    ic_type = 'step'
    solver.ic_type = ic_type  # 保存初始条件类型用于精确解计算
    u_final, history = solver.solve(T, boundary='periodic', ic_type=ic_type, verbose=True)
    
    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 子图1: 初始条件和最终解对比
    ax1 = axes[0, 0]
    ax1.plot(solver.x, solver.u0, 'b-', linewidth=2, label='初始条件')
    ax1.plot(solver.x, u_final, 'r-', linewidth=2, label=f'数值解 (t={T})')
    
    # 精确解
    exact_final = solver.exact_solution(T, ic_type)
    ax1.plot(solver.x, exact_final, 'g--', linewidth=2, label='精确解', alpha=0.7)
    
    ax1.set_xlabel('位置 x')
    ax1.set_ylabel('u(x,t)')
    ax1.set_title(f'Crank-Nicolson格式 (CFL={solver.cfl:.2f})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 子图2: 误差分布
    ax2 = axes[0, 1]
    error = u_final - exact_final
    ax2.plot(solver.x, error, 'r-', linewidth=2)
    ax2.fill_between(solver.x, 0, error, alpha=0.3, color='red')
    ax2.set_xlabel('位置 x')
    ax2.set_ylabel('误差')
    ax2.set_title(f'误差分布 (t={T})')
    ax2.grid(True, alpha=0.3)
    
    # 计算误差统计
    error_L1, error_L2, error_Linf = solver.compute_error(T, u_final)
    error_text = f'L1误差: {error_L1:.2e}\nL2误差: {error_L2:.2e}\n最大误差: {error_Linf:.2e}'
    ax2.text(0.02, 0.98, error_text, transform=ax2.transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 子图3: 时间演化（动画的几帧）
    ax3 = axes[1, 0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(history)))
    
    for i, (t, u) in enumerate(history):
        alpha = 0.3 + 0.7 * (i / len(history))
        ax3.plot(solver.x, u, '-', color=colors[i], alpha=alpha, 
                label=f't={t:.2f}' if i % 3 == 0 else "")
    
    ax3.set_xlabel('位置 x')
    ax3.set_ylabel('u(x,t)')
    ax3.set_title('时间演化')
    ax3.legend(loc='upper right', fontsize=8)
    ax3.grid(True, alpha=0.3)
    
    # 子图4: 不同CFL数对比
    ax4 = axes[1, 1]
    cfl_list = [0.1, 0.5, 1.0, 2.0, 5.0]
    colors_cfl = plt.cm.plasma(np.linspace(0, 1, len(cfl_list)))
    
    for idx, cfl in enumerate(cfl_list):
        # 创建新求解器
        solver_test = CrankNicolson1D(a=1.0, L=1.0, nx=101, cfl=cfl)
        solver_test.set_initial_condition(ic_type)
        
        # 只计算一个时间步，看耗散/色散
        u_test, _ = solver_test.solve(0.1, boundary='periodic', ic_type=ic_type, verbose=False)
        exact_test = solver_test.exact_solution(0.1, ic_type)
        
        ax4.plot(solver_test.x, u_test - exact_test, '-', color=colors_cfl[idx], 
                linewidth=1.5, label=f'CFL={cfl}')
    
    ax4.set_xlabel('位置 x')
    ax4.set_ylabel('误差')
    ax4.set_title('不同CFL数下的误差对比 (t=0.1)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(f'crank_nicolson_step_cfl{solver.cfl:.1f}.png', dpi=150, bbox_inches='tight')
    
    plt.show()
    
    # 打印误差统计
    print("\n" + "="*50)
    print("误差分析:")
    print(f"L1误差 (积分绝对误差): {error_L1:.6e}")
    print(f"L2误差 (均方根误差): {error_L2:.6e}")
    print(f"L∞误差 (最大绝对误差): {error_Linf:.6e}")
    print("="*50)
    
    return error_L1, error_L2, error_Linf


def convergence_study():
    """
    网格收敛性研究
    """
    print("进行网格收敛性研究...")
    
    # 测试不同的网格分辨率
    nx_list = [51, 101, 201, 401, 801]
    errors_L1 = []
    errors_L2 = []
    errors_Linf = []
    
    T = 0.5  # 最终时间
    cfl = 0.5  # 固定CFL数
    
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    for nx in nx_list:
        print(f"\n网格数 nx = {nx}")
        
        # 创建求解器
        solver = CrankNicolson1D(a=1.0, L=1.0, nx=nx, cfl=cfl)
        solver.ic_type = 'step'
        
        # 求解
        u_final, _ = solver.solve(T, boundary='periodic', ic_type='step', verbose=False)
        
        # 计算误差
        error_L1, error_L2, error_Linf = solver.compute_error(T, u_final)
        errors_L1.append(error_L1)
        errors_L2.append(error_L2)
        errors_Linf.append(error_Linf)
        
        print(f"  误差: L1={error_L1:.2e}, L2={error_L2:.2e}, L∞={error_Linf:.2e}")
        
        # 绘制该分辨率的解
        if nx in [101, 401]:
            ax[0].plot(solver.x, u_final, '-', linewidth=1.5, 
                      label=f'nx={nx}', alpha=0.7)
    
    # 绘制精确解
    solver_ref = CrankNicolson1D(a=1.0, L=1.0, nx=801, cfl=cfl)
    exact = solver_ref.exact_solution(T, 'step')
    ax[0].plot(solver_ref.x, exact, 'k--', linewidth=2, label='精确解', alpha=0.5)
    
    ax[0].set_xlabel('位置 x')
    ax[0].set_ylabel('u(x,t)')
    ax[0].set_title(f'不同网格分辨率的解 (t={T})')
    ax[0].legend()
    ax[0].grid(True, alpha=0.3)
    
    # 绘制收敛率
    dx_list = [1.0/(nx-1) for nx in nx_list]
    
    # 拟合收敛率
    log_dx = np.log(np.array(dx_list))
    log_L1 = np.log(np.array(errors_L1))
    log_L2 = np.log(np.array(errors_L2))
    log_Linf = np.log(np.array(errors_Linf))
    
    # 最小二乘拟合斜率（收敛阶）
    p_L1 = np.polyfit(log_dx, log_L1, 1)[0]
    p_L2 = np.polyfit(log_dx, log_L2, 1)[0]
    p_Linf = np.polyfit(log_dx, log_Linf, 1)[0]
    
    ax[1].loglog(dx_list, errors_L1, 'o-', linewidth=2, markersize=8, label=f'L1误差 (阶数={p_L1:.2f})')
    ax[1].loglog(dx_list, errors_L2, 's-', linewidth=2, markersize=8, label=f'L2误差 (阶数={p_L2:.2f})')
    ax[1].loglog(dx_list, errors_Linf, '^-', linewidth=2, markersize=8, label=f'L∞误差 (阶数={p_Linf:.2f})')
    
    # 绘制参考线（一阶和二阶收敛）
    ref_dx = np.array([dx_list[0], dx_list[-1]])
    ax[1].loglog(ref_dx, 0.1*ref_dx, 'k--', linewidth=1, alpha=0.5, label='一阶收敛')
    ax[1].loglog(ref_dx, 0.01*ref_dx**2, 'k:', linewidth=1, alpha=0.5, label='二阶收敛')
    
    ax[1].set_xlabel('空间步长 Δx')
    ax[1].set_ylabel('误差')
    ax[1].set_title('收敛性分析')
    ax[1].legend()
    ax[1].grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('convergence_study.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("\n" + "="*50)
    print("收敛性分析结果:")
    print(f"L1误差收敛阶: {p_L1:.3f}")
    print(f"L2误差收敛阶: {p_L2:.3f}")
    print(f"L∞误差收敛阶: {p_Linf:.3f}")
    print("="*50)
    
    return p_L1, p_L2, p_Linf


def compare_with_explicit():
    """
    对比Crank-Nicolson和显式格式
    """
    print("对比Crank-Nicolson和显式格式...")
    
    # 参数
    nx = 101
    L = 1.0
    a = 1.0
    T = 0.5
    
    # 创建求解器
    cn_solver = CrankNicolson1D(a=a, L=L, nx=nx, cfl=1.0)  # CFL=1.0
    cn_solver.set_initial_condition('step')
    
    # 显式欧拉求解器
    class ExplicitEuler1D:
        def __init__(self, a=1.0, L=1.0, nx=101, cfl=0.5):
            self.a = a
            self.L = L
            self.nx = nx
            self.dx = L / (nx - 1)
            self.x = np.linspace(0, L, nx)
            self.dt = cfl * self.dx / abs(a)
            self.cfl = cfl
        
        def set_initial_condition(self, ic_type='step'):
            if ic_type == 'step':
                self.u = np.where(self.x < 0.5, 1.0, 0.0)
            self.u0 = self.u.copy()
        
        def step(self, u_old, boundary='periodic'):
            u_new = u_old.copy()
            n = len(u_old)
            
            if boundary == 'periodic':
                for i in range(n):
                    left = u_old[(i-1) % n]
                    right = u_old[(i+1) % n]
                    # 中心差分 + 显式欧拉
                    u_new[i] = u_old[i] - self.a * self.dt / (2 * self.dx) * (right - left)
            
            return u_new
        
        def solve(self, T, ic_type='step'):
            self.set_initial_condition(ic_type)
            u = self.u.copy()
            
            nt = int(np.ceil(T / self.dt))
            for n in range(nt):
                u = self.step(u, 'periodic')
            
            return u
    
    # 创建显式求解器（CFL必须很小！）
    ee_solver = ExplicitEuler1D(a=a, L=L, nx=nx, cfl=0.1)  # 小CFL
    
    # 求解
    u_cn, _ = cn_solver.solve(T, boundary='periodic', ic_type='step', verbose=False)
    u_ee = ee_solver.solve(T, ic_type='step')
    
    # 精确解
    exact = cn_solver.exact_solution(T, 'step')
    
    # 绘制对比
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 子图1: 解对比
    ax1 = axes[0, 0]
    ax1.plot(cn_solver.x, cn_solver.u0, 'k-', linewidth=2, label='初始条件')
    ax1.plot(cn_solver.x, u_cn, 'b-', linewidth=2, label=f'Crank-Nicolson (CFL={cn_solver.cfl:.1f})')
    ax1.plot(cn_solver.x, u_ee, 'r--', linewidth=2, label=f'显式欧拉 (CFL={ee_solver.cfl:.1f})')
    ax1.plot(cn_solver.x, exact, 'g:', linewidth=3, label='精确解', alpha=0.7)
    
    ax1.set_xlabel('位置 x')
    ax1.set_ylabel('u(x,t)')
    ax1.set_title(f'方法对比 (t={T})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 子图2: 误差对比
    ax2 = axes[0, 1]
    error_cn = u_cn - exact
    error_ee = u_ee - exact
    
    ax2.plot(cn_solver.x, error_cn, 'b-', linewidth=2, label='Crank-Nicolson误差')
    ax2.plot(cn_solver.x, error_ee, 'r--', linewidth=2, label='显式欧拉误差')
    
    ax2.set_xlabel('位置 x')
    ax2.set_ylabel('误差')
    ax2.set_title('误差分布')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 子图3: 不同CFL下的Crank-Nicolson
    ax3 = axes[1, 0]
    cfl_list_cn = [0.5, 1.0, 2.0, 5.0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(cfl_list_cn)))
    
    for idx, cfl in enumerate(cfl_list_cn):
        solver = CrankNicolson1D(a=a, L=L, nx=nx, cfl=cfl)
        solver.set_initial_condition('step')
        u, _ = solver.solve(T, boundary='periodic', ic_type='step', verbose=False)
        
        error = u - solver.exact_solution(T, 'step')
        ax3.plot(solver.x, error, '-', color=colors[idx], linewidth=1.5, label=f'CFL={cfl}')
    
    ax3.set_xlabel('位置 x')
    ax3.set_ylabel('误差')
    ax3.set_title('Crank-Nicolson不同CFL数的误差')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 子图4: 不同CFL下的显式欧拉（如果稳定的话）
    ax4 = axes[1, 1]
    cfl_list_ee = [0.1, 0.2, 0.3, 0.4]
    colors_ee = plt.cm.plasma(np.linspace(0, 1, len(cfl_list_ee)))
    
    for idx, cfl in enumerate(cfl_list_ee):
        try:
            solver = ExplicitEuler1D(a=a, L=L, nx=nx, cfl=cfl)
            solver.set_initial_condition('step')
            u = solver.solve(T, ic_type='step')
            
            error = u - exact
            ax4.plot(solver.x, error, '-', color=colors_ee[idx], linewidth=1.5, label=f'CFL={cfl}')
        except Exception as e:
            print(f"CFL={cfl}时显式欧拉失败: {e}")
    
    ax4.set_xlabel('位置 x')
    ax4.set_ylabel('误差')
    ax4.set_title('显式欧拉不同CFL数的误差')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('comparison_cn_vs_explicit.png', dpi=150, bbox_inches='tight')
    plt.show()


# 主程序
if __name__ == "__main__":
    print("="*60)
    print("一维线性对流方程的Crank-Nicolson格式求解")
    print("方程: ∂u/∂t + a ∂u/∂x = 0")
    print("初始条件: 阶跃函数")
    print("="*60)
    
    # 示例1: 基本测试
    print("\n1. 基本测试 (CFL=0.5)...")
    solver = CrankNicolson1D(a=1.0, L=1.0, nx=201, cfl=0.5)
    plot_results(solver, T=0.5, save_fig=True)
    
    # 示例2: 收敛性研究
    print("\n2. 网格收敛性研究...")
    convergence_study()
    
    # 示例3: 与显式方法对比
    print("\n3. 与显式欧拉方法对比...")
    compare_with_explicit()
    
    print("\n" + "="*60)
    print("所有测试完成!")
    print("="*60)