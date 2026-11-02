import numpy as np
import matplotlib.pyplot as plt
import re
import os

def get_order_from_filename(filename):
    """从文件名提取阶数，默认返回3"""
    match = re.search(r'_order(\d+)', filename)
    return int(match.group(1)) if match else 3

def read_results(filename):
    """读取结果文件"""
    try:
        data = np.loadtxt(filename, comments='#')
        return data[:, 0], data[:, 1]
    except:
        return np.array([]), np.array([])

def find_files():
    """查找结果文件，优先使用带order标记的文件"""
    files = [f for f in os.listdir('.') if f.endswith('.txt')]
    
    # 查找带order标记的文件
    eno_files = [f for f in files if f.startswith('eno_results')]
    weno_files = [f for f in files if f.startswith('weno_results')]
    analytical_files = [f for f in files if f.startswith('analytical_results')]
    
    # 选择文件
    def select_file(file_list):
        # 优先选择带order标记的文件
        order_files = [f for f in file_list if '_order' in f]
        if order_files:
            return order_files[0]
        elif 'results.txt' in file_list:  # 默认文件
            return 'results.txt'
        return None
    
    eno_file = select_file(eno_files)
    weno_file = select_file(weno_files)
    analytical_file = select_file(analytical_files)
    
    if not all([eno_file, weno_file, analytical_file]):
        # 尝试默认文件名
        defaults = ['eno_results.txt', 'weno_results.txt', 'analytical_results.txt']
        return defaults[0], defaults[1], defaults[2], 3
    
    # 从ENO文件名提取阶数
    order = get_order_from_filename(eno_file)
    
    return eno_file, weno_file, analytical_file, order

def main():
    print("Post-processing Fortran output...")
    
    # 查找文件并确定阶数
    eno_file, weno_file, analytical_file, order = find_files()
    
    print(f"Files: {eno_file}, {weno_file}, {analytical_file}")
    print(f"Order: {order}")
    
    # 读取数据
    x_eno, u_eno = read_results(eno_file)
    x_weno, u_weno = read_results(weno_file)
    x_analytical, u_analytical = read_results(analytical_file)
    
    if len(x_eno) == 0:
        print("Error: No data found.")
        return
    
    # 绘图
    plt.figure(figsize=(10, 6))
    
    # 动态标签
    plt.plot(x_eno, u_eno, 'bo-', linewidth=1, markersize=5, 
             markerfacecolor='none', label=f'ENO{order}')
    plt.plot(x_weno, u_weno, 'gs-', linewidth=1, markersize=5,
             markerfacecolor='none', label=f'WENO{order}')
    plt.plot(x_analytical, u_analytical, 'r-', linewidth=2, label='Analytical')
    
    # 动态标题
    plt.title(f'1D Convection: ENO{order} vs WENO{order}')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 保存文件
    output_file = f'comparison_order{order}.png'
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved as {output_file}")
    plt.show()

if __name__ == "__main__":
    main()