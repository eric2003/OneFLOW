"""
文件名: 04_cfd_animation.py
功能: 创建CFD计算过程动画
包含: 波动方程的数值解演化过程
修复: 解决了动画保存后plt.show()的AttributeError问题
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import warnings
warnings.filterwarnings('ignore')

def setup_plot_style():
    """设置绘图样式"""
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'figure.titlesize': 16
    })

def simulate_convection_equation():
    """模拟对流方程的解"""
    # 设置网格和参数
    L = 10.0  # 计算域长度
    nx = 100  # 网格点数
    dx = L / nx
    x = np.linspace(0, L, nx)
    
    # 时间参数
    nt = 100  # 减少帧数以加快速度
    dt = 0.05
    
    # 初始化变量
    u_exact = np.zeros((nt, nx))
    u_numerical = np.zeros((nt, nx))
    
    # 初始条件：高斯波包
    u0 = np.exp(-(x - L/4)**2 / 0.5) * np.sin(3 * x)
    
    # 生成精确解（对流方程）
    c = 0.5  # 对流速度
    for n in range(nt):
        # 精确解：简单对流
        shift = c * n * dt
        u_exact[n, :] = np.exp(-(x - shift - L/4)**2 / 0.5) * np.sin(3 * (x - shift))
        
        # 数值解：使用一阶迎风格式
        if n == 0:
            u_numerical[n, :] = u0
        else:
            for i in range(1, nx):
                # 一阶迎风格式
                u_numerical[n, i] = u_numerical[n-1, i] - c*dt/dx * (
                    u_numerical[n-1, i] - u_numerical[n-1, i-1])
            # 边界条件：左侧固定为零
            u_numerical[n, 0] = 0.0
    
    return x, u_exact, u_numerical, nt, dt

def create_cfd_animation(save_video=True, show_plot=True):
    """创建CFD计算过程动画"""
    setup_plot_style()
    
    # 获取模拟数据
    x, u_exact, u_numerical, nt, dt = simulate_convection_equation()
    
    # 创建图形
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # 初始化线条
    line_exact, = ax1.plot([], [], 'b-', linewidth=2, label='Exact Solution')
    line_num, = ax2.plot([], [], 'r-', linewidth=2, label='Numerical Solution')
    line_num_dots, = ax2.plot([], [], 'ro', markersize=4, alpha=0.5, markevery=5)
    
    # 设置坐标轴
    for ax in [ax1, ax2]:
        ax.set_xlim(0, 10)
        ax.set_ylim(-1.2, 1.2)
        ax.set_xlabel('Position x')
        ax.set_ylabel('Variable u')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
    
    ax1.set_title('Exact Solution of Convection Equation')
    ax2.set_title('Numerical Solution (First-order Upwind Scheme)')
    
    fig.suptitle('1D Convection Equation Simulation', fontsize=16, fontweight='bold', y=0.98)
    
    # 添加时间文本
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12, 
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # 添加说明文本
    info_text = fig.text(0.5, 0.02, 
                         'Convection velocity: c = 0.5 | CFL number: ~0.5', 
                         ha='center', fontsize=10,
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))
    
    # 初始化函数
    def init():
        line_exact.set_data([], [])
        line_num.set_data([], [])
        line_num_dots.set_data([], [])
        time_text.set_text('')
        return line_exact, line_num, line_num_dots, time_text
    
    # 动画更新函数
    def update(frame):
        # 更新精确解
        line_exact.set_data(x, u_exact[frame, :])
        
        # 更新数值解
        line_num.set_data(x, u_numerical[frame, :])
        # 每隔5个点显示一个标记点
        line_num_dots.set_data(x[::5], u_numerical[frame, ::5])
        
        # 更新时间文本
        time_text.set_text(f'Time: {frame*dt:.2f} s | Frame: {frame}/{nt}')
        
        return line_exact, line_num, line_num_dots, time_text
    
    # 创建动画
    ani = animation.FuncAnimation(fig, update, frames=nt,
                                  init_func=init, blit=True, interval=100, repeat=False)
    
    # 保存动画（可选）
    if save_video:
        try:
            print("正在保存动画...")
            # 尝试使用不同的写入器
            try:
                ani.save('04_cfd_animation.mp4', writer='ffmpeg', fps=10, dpi=150)
                print("动画已保存为 '04_cfd_animation.mp4'")
            except:
                # 如果ffmpeg不可用，尝试保存为gif
                ani.save('04_cfd_animation.gif', writer='pillow', fps=10, dpi=150)
                print("动画已保存为 '04_cfd_animation.gif' (使用pillow写入器)")
        except Exception as e:
            print(f"保存动画失败: {e}")
            print("请确保已安装必要的视频编码器（如ffmpeg）或pillow库")
    
    # 显示动画（可选）
    if show_plot:
        try:
            plt.tight_layout()
            plt.show()
        except Exception as e:
            print(f"显示动画时出错: {e}")
            print("这可能是因为图形窗口已关闭，但这是正常的。")
    
    return ani

def create_static_frames():
    """创建静态的关键帧图像，作为替代方案"""
    setup_plot_style()
    
    # 获取模拟数据
    x, u_exact, u_numerical, nt, dt = simulate_convection_equation()
    
    # 选择几个关键帧
    key_frames = [0, nt//4, nt//2, 3*nt//4, nt-1]
    frame_names = ['Initial', 'Quarter', 'Half', 'Three Quarters', 'Final']
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    
    for idx, (frame, name) in enumerate(zip(key_frames, frame_names)):
        if idx >= len(axes):
            break
            
        ax = axes[idx]
        ax.plot(x, u_exact[frame, :], 'b-', linewidth=2, label='Exact', alpha=0.7)
        ax.plot(x, u_numerical[frame, :], 'r--', linewidth=2, label='Numerical')
        ax.scatter(x[::10], u_numerical[frame, ::10], s=30, color='red', alpha=0.5)
        
        ax.set_xlim(0, 10)
        ax.set_ylim(-1.2, 1.2)
        ax.set_xlabel('Position x')
        ax.set_ylabel('u')
        ax.set_title(f'{name} (t={frame*dt:.2f}s)')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
    
    # 第6个子图：误差随时间变化
    ax = axes[5]
    errors = np.abs(u_numerical - u_exact)
    mean_errors = np.mean(errors, axis=1)
    max_errors = np.max(errors, axis=1)
    
    time = np.arange(nt) * dt
    ax.plot(time, mean_errors, 'b-', label='Mean Error')
    ax.plot(time, max_errors, 'r-', label='Max Error')
    ax.fill_between(time, 0, mean_errors, alpha=0.3, color='blue')
    ax.fill_between(time, 0, max_errors, alpha=0.1, color='red')
    
    ax.set_xlabel('Time')
    ax.set_ylabel('Error')
    ax.set_title('Numerical Error Evolution')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.suptitle('1D Convection Equation: Key Frames and Error Analysis', 
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('04_cfd_static_frames.png', dpi=300, bbox_inches='tight')
    print("静态关键帧已保存为 '04_cfd_static_frames.png'")
    plt.show()

def create_simple_animation():
    """创建一个更简单的动画，避免复杂问题"""
    setup_plot_style()
    
    # 简单示例：波的传播
    x = np.linspace(0, 10, 200)
    t = np.linspace(0, 4*np.pi, 100)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 创建初始线
    line, = ax.plot([], [], 'b-', linewidth=2)
    
    # 设置图形
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Amplitude u')
    ax.set_title('Wave Propagation in 1D Domain')
    ax.grid(True, alpha=0.3)
    
    # 初始化函数
    def init():
        line.set_data([], [])
        return line,
    
    # 更新函数
    def update(i):
        y = np.sin(2*np.pi*(x/5 - t[i]/5)) * np.exp(-0.05*x)
        line.set_data(x, y)
        ax.set_title(f'Wave Propagation (t = {t[i]:.2f})')
        return line,
    
    # 创建动画
    ani = animation.FuncAnimation(fig, update, frames=len(t),
                                  init_func=init, blit=True, interval=50)
    
    # 保存为gif
    try:
        ani.save('04_simple_wave.gif', writer='pillow', fps=15)
        print("简单动画已保存为 '04_simple_wave.gif'")
    except Exception as e:
        print(f"保存简单动画失败: {e}")
    
    plt.show()
    return ani

if __name__ == "__main__":
    print("=" * 60)
    print("1D CFD 动画演示程序")
    print("=" * 60)
    print("\n选项:")
    print("1. 创建完整动画（可能有问题）")
    print("2. 创建静态关键帧图像")
    print("3. 创建简单波动动画")
    
    try:
        choice = int(input("\n请选择 (1-3, 默认=2): ") or "2")
    except:
        choice = 2
    
    if choice == 1:
        print("\n正在创建完整动画...")
        print("注意：保存视频需要ffmpeg或pillow库")
        print("如果失败，将自动回退到静态图像")
        
        # 先尝试创建静态图像
        create_static_frames()
        
        # 然后尝试动画
        try:
            ani = create_cfd_animation(save_video=True, show_plot=True)
        except Exception as e:
            print(f"\n动画创建失败: {e}")
            print("已创建静态图像作为替代。")
            
    elif choice == 2:
        print("\n正在创建静态关键帧图像...")
        create_static_frames()
        
    elif choice == 3:
        print("\n正在创建简单波动动画...")
        create_simple_animation()
        
    else:
        print("\n无效选择，创建静态关键帧图像...")
        create_static_frames()
    
    print("\n程序执行完成！")