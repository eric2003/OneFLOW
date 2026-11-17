"""
文件名: 04_cfd_animation_fixed.py
功能: 完全按照简单模式修复的CFD动画
参考: 你的成功示例代码
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def create_simple_animation():
    """创建简单的CFD动画（只显示）"""
    print("正在创建CFD动画...")
    
    # 创建图形和轴
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Velocity u')
    ax.grid(True, alpha=0.3)
    
    # 初始化数据
    x = np.linspace(0, 10, 200)
    line_exact, = ax.plot(x, np.zeros_like(x), 'b-', linewidth=2, label='Exact Solution')
    line_num, = ax.plot(x, np.zeros_like(x), 'r--', linewidth=2, label='Numerical Solution')
    ax.legend()
    
    # 添加标题
    ax.set_title('1D Convection Equation Simulation')
    
    # 动画更新函数
    def animate(frame):
        t = frame * 0.05  # 时间步长
        
        # 精确解：波动传播
        u_exact = np.sin(2 * np.pi * (x/5 - t/2)) * np.exp(-0.1*x)
        
        # 数值解：添加一些数值误差（模拟迎风格式的耗散）
        u_num = u_exact * (1 - 0.05*t) + 0.1*np.sin(5*x)*np.sin(t)
        
        line_exact.set_ydata(u_exact)
        line_num.set_ydata(u_num)
        ax.set_title(f'1D Convection Equation - Time: {t:.2f}s')
        
        return line_exact, line_num,
    
    # 创建动画：100帧，每帧间隔50ms
    ani = animation.FuncAnimation(fig, animate, frames=100, interval=50, blit=True)
    
    print("动画准备就绪，显示窗口...")
    print("提示：关闭窗口以退出程序。")
    
    plt.tight_layout()
    plt.show()
    
    return ani

def create_and_save_animation():
    """创建并保存动画（先保存再显示分开处理）"""
    print("正在创建并保存CFD动画...")
    
    # ========== 第一部分：只保存动画 ==========
    print("\n1. 保存动画到文件...")
    
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set_xlim(0, 10)
    ax1.set_ylim(-1.5, 1.5)
    
    x = np.linspace(0, 10, 200)
    line1, = ax1.plot(x, np.zeros_like(x), 'b-', lw=2)
    line2, = ax1.plot(x, np.zeros_like(x), 'r--', lw=2)
    
    def update_save(frame):
        t = frame * 0.05
        u_exact = np.sin(2 * np.pi * (x/5 - t/2))
        u_num = u_exact * (1 - 0.05*t)
        
        line1.set_ydata(u_exact)
        line2.set_ydata(u_num)
        ax1.set_title(f'Time: {t:.2f}s')
        
        return line1, line2,
    
    ani_save = animation.FuncAnimation(fig1, update_save, frames=50, blit=True)
    
    # 保存动画
    try:
        ani_save.save('04_cfd_save_only.gif', writer='pillow', fps=15, dpi=100)
        print("✓ 动画已保存为 '04_cfd_save_only.gif'")
    except Exception as e:
        print(f"保存失败: {e}")
    
    # 重要：关闭保存用的图形
    plt.close(fig1)
    
    # ========== 第二部分：只显示动画 ==========
    print("\n2. 显示动画...")
    create_simple_animation()

def create_static_images():
    """创建静态图像"""
    print("正在创建静态图像...")
    
    x = np.linspace(0, 10, 200)
    times = [0, 0.5, 1.0, 1.5, 2.0]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()
    
    for i, (ax, t) in enumerate(zip(axes[:5], times)):
        u_exact = np.sin(2 * np.pi * (x/5 - t/2))
        u_num = u_exact * (1 - 0.05*t)
        
        ax.plot(x, u_exact, 'b-', lw=2, label='Exact')
        ax.plot(x, u_num, 'r--', lw=2, label='Numerical')
        ax.set_xlim(0, 10)
        ax.set_ylim(-1.5, 1.5)
        ax.set_xlabel('Position x')
        ax.set_ylabel('Velocity u')
        ax.set_title(f'Time = {t:.1f}s')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # 最后一个子图显示误差
    ax = axes[5]
    t_vals = np.linspace(0, 2, 50)
    errors = []
    
    for t in t_vals:
        u_exact = np.sin(2 * np.pi * (x/5 - t/2))
        u_num = u_exact * (1 - 0.05*t)
        error = np.max(np.abs(u_exact - u_num))
        errors.append(error)
    
    ax.plot(t_vals, errors, 'k-', lw=2)
    ax.fill_between(t_vals, 0, errors, alpha=0.3)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Max Error')
    ax.set_title('Numerical Error over Time')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('1D Convection Equation - Static Frames', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig('04_cfd_static.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✓ 静态图像已保存为 '04_cfd_static.png'")

def main_menu():
    """主菜单"""
    print("=" * 60)
    print("CFD 动画演示程序")
    print("=" * 60)
    print("\n选项:")
    print("1. 显示动画（简单，不会出错）")
    print("2. 保存并显示动画（分开处理）")
    print("3. 创建静态图像")
    
    try:
        choice = input("\n请选择 (1-3, 默认=1): ").strip()
        if choice == "":
            choice = "1"
        choice = int(choice)
    except:
        choice = 1
    
    if choice == 1:
        create_simple_animation()
    elif choice == 2:
        create_and_save_animation()
    elif choice == 3:
        create_static_images()
    else:
        create_simple_animation()

if __name__ == "__main__":
    main_menu()