"""
文件名: 04_cfd_animation.py
功能: 创建CFD计算过程动画
修复: 分离保存动画和显示动画的过程，避免AttributeError
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import warnings
import sys
import os

warnings.filterwarnings('ignore')

def setup_plot_style():
    """设置绘图样式"""
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'figure.dpi': 100,
        'savefig.dpi': 150,
        'animation.embed_limit': 100  # 增加嵌入限制
    })

def simulate_convection_equation():
    """模拟对流方程的解"""
    # 设置网格和参数
    L = 10.0  # 计算域长度
    nx = 100  # 网格点数
    dx = L / nx
    x = np.linspace(0, L, nx)
    
    # 时间参数
    nt = 80  # 减少帧数以加快速度
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

def save_animation_only():
    """仅保存动画，不显示"""
    print("正在创建并保存动画...")
    
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
        line_num_dots.set_data(x[::5], u_numerical[frame, ::5])
        
        # 更新时间文本
        time_text.set_text(f'Time: {frame*dt:.2f} s')
        
        return line_exact, line_num, line_num_dots, time_text
    
    # 创建动画对象
    ani = FuncAnimation(fig, update, frames=nt,
                       init_func=init, blit=True)
    
    # 尝试保存动画
    saved = False
    
    # 首先尝试保存为mp4
    try:
        print("尝试保存为MP4格式...")
        writer = FFMpegWriter(fps=10, metadata=dict(artist='CFD Visualization'), bitrate=1800)
        ani.save('04_cfd_animation.mp4', writer=writer)
        print("✓ 动画已成功保存为 '04_cfd_animation.mp4'")
        saved = True
    except Exception as e:
        print(f"MP4保存失败: {e}")
    
    # 如果mp4失败，尝试保存为gif
    if not saved:
        try:
            print("尝试保存为GIF格式...")
            writer = PillowWriter(fps=10)
            ani.save('04_cfd_animation.gif', writer=writer)
            print("✓ 动画已保存为 '04_cfd_animation.gif'")
            saved = True
        except Exception as e:
            print(f"GIF保存失败: {e}")
    
    # 关闭图形以释放资源
    plt.close(fig)
    
    if saved:
        print("\n动画保存完成！文件已生成。")
        print("您可以使用视频播放器查看保存的动画文件。")
    else:
        print("\n无法保存动画，将创建静态图像...")
        create_static_frames()
    
    return saved

def create_static_frames():
    """创建静态的关键帧图像"""
    print("正在创建静态关键帧图像...")
    
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
    print("✓ 静态关键帧已保存为 '04_cfd_static_frames.png'")
    plt.close(fig)
    
    return True

def show_animation_only():
    """仅显示动画，不保存"""
    print("正在创建动画以供显示...")
    
    # 获取模拟数据
    x, u_exact, u_numerical, nt, dt = simulate_convection_equation()
    
    # 创建新的图形
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
    
    fig.suptitle('1D Convection Equation Simulation (Live)', fontsize=16, fontweight='bold', y=0.98)
    
    # 添加时间文本
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=12, 
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
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
        line_num_dots.set_data(x[::5], u_numerical[frame, ::5])
        
        # 更新时间文本
        time_text.set_text(f'Time: {frame*dt:.2f} s | Frame: {frame}/{nt}')
        
        return line_exact, line_num, line_num_dots, time_text
    
    # 创建动画
    ani = FuncAnimation(fig, update, frames=nt,
                       init_func=init, blit=True, interval=100, repeat=False)
    
    print("动画准备就绪，显示窗口...")
    print("提示：关闭窗口以继续程序。")
    
    plt.tight_layout()
    plt.show()
    
    return ani

def check_ffmpeg_available():
    """检查ffmpeg是否可用"""
    try:
        import subprocess
        result = subprocess.run(['ffmpeg', '-version'], 
                               capture_output=True, text=True, shell=True)
        return result.returncode == 0
    except:
        return False

def main():
    """主函数"""
    print("=" * 60)
    print("1D CFD 动画演示程序")
    print("=" * 60)
    
    # 检查ffmpeg
    has_ffmpeg = check_ffmpeg_available()
    print(f"FFmpeg可用: {'✓' if has_ffmpeg else '✗'}")
    
    print("\n选项:")
    print("1. 保存动画为视频文件（不显示）")
    print("2. 显示动画（不保存）")
    print("3. 创建静态关键帧图像")
    
    try:
        choice = int(input("\n请选择 (1-3, 默认=1): ") or "1")
    except:
        choice = 1
    
    setup_plot_style()
    
    if choice == 1:
        print("\n" + "=" * 50)
        print("选项1: 保存动画为视频文件")
        print("=" * 50)
        save_animation_only()
        
    elif choice == 2:
        print("\n" + "=" * 50)
        print("选项2: 显示动画")
        print("=" * 50)
        print("注意：动画将在新窗口中显示")
        print("关闭窗口后程序将继续运行")
        show_animation_only()
        
    elif choice == 3:
        print("\n" + "=" * 50)
        print("选项3: 创建静态关键帧图像")
        print("=" * 50)
        create_static_frames()
        
    else:
        print("\n无效选择，使用默认选项...")
        save_animation_only()
    
    print("\n" + "=" * 60)
    print("程序执行完成！")
    print("生成的文件:")
    
    # 列出生成的文件
    files_to_check = [
        ('04_cfd_animation.mp4', '动画视频文件'),
        ('04_cfd_animation.gif', '动画GIF文件'),
        ('04_cfd_static_frames.png', '静态关键帧图像')
    ]
    
    for filename, description in files_to_check:
        if os.path.exists(filename):
            file_size = os.path.getsize(filename) / 1024  # KB
            print(f"  ✓ {filename} - {description} ({file_size:.1f} KB)")
        else:
            print(f"  ✗ {filename} - 未生成")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n程序被用户中断。")
    except Exception as e:
        print(f"\n程序执行出错: {e}")
        import traceback
        traceback.print_exc()