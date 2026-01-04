#!/usr/bin/env python3
# python/plot_results.py
"""
Fortran CFD 结果可视化脚本
与Julia的plotter.jl功能类似，但直接读取Fortran生成的文本文件
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob
from pathlib import Path
import re

class FortranResultsPlotter:
    """Fortran结果绘图器"""
    
    def __init__(self, style='default'):
        self.style = style
        self.setup_styles()
    
    def setup_styles(self):
        """设置绘图样式"""
        if self.style == 'default':
            self.styles = {
                'numerical': {
                    'color': 'blue',
                    'linestyle': '-',
                    'marker': 'o',
                    'markerfacecolor': 'none',
                    'markersize': 4,
                    'linewidth': 1
                },
                'analytical': {
                    'color': 'red',
                    'linestyle': '--',
                    'marker': '',
                    'linewidth': 1.5
                },
                'comparison': [
                    {'color': 'black', 'linestyle': '-', 'marker': 'o', 'markerfacecolor': 'none'},
                    {'color': 'blue', 'linestyle': '--', 'marker': 's', 'markerfacecolor': 'none'},
                    {'color': 'green', 'linestyle': ':', 'marker': '^', 'markerfacecolor': 'none'}
                ]
            }
    
    def read_fortran_results(self, filename):
        """读取Fortran生成的结果文件"""
        data = {}
        
        try:
            with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
            
            data_lines = []
            in_data_section = False  # 新状态：是否在数据区域
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                
                # 跳过所有分隔线
                if line.startswith("===="):
                    continue
                
                # 检测数据区域开始
                if line == "DATA: x, numerical, analytical":
                    in_data_section = True
                    continue
                
                if not in_data_section:
                    # 解析头部信息
                    if "Solver:" in line:
                        data['solver'] = line.split(":", 1)[1].strip()
                    elif "Scheme:" in line:
                        data['scheme'] = line.split(":", 1)[1].strip()
                    elif "Order:" in line:
                        if "RK Order:" in line:
                            data['rk_order'] = int(line.split(":", 1)[1].strip())
                        else:
                            data['order'] = int(line.split(":", 1)[1].strip())
                    elif "Current Time:" in line:
                        data['time'] = float(line.split(":", 1)[1].strip())
                    elif "Grid Points:" in line:
                        data['n_points'] = int(line.split(":", 1)[1].strip())
                else:
                    # 解析数据行
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            x = float(parts[0])
                            numerical = float(parts[1])
                            analytical = float(parts[2])
                            data_lines.append([x, numerical, analytical])
                        except ValueError:
                            continue  # 忽略无法解析的行
            
            if data_lines:
                import numpy as np
                data_array = np.array(data_lines)
                data['x'] = data_array[:, 0]
                data['numerical'] = data_array[:, 1]
                data['analytical'] = data_array[:, 2]
                
                print(f"Read {len(data['x'])} points from {filename}")
                print(f"  Solver: {data.get('solver', 'N/A')}")
                print(f"  Scheme: {data.get('scheme', 'N/A')} order {data.get('order', 'N/A')}")
                print(f"  Time: {data.get('time', 'N/A')}")
                
                return data
            else:
                print(f"Warning: No data found in {filename}")
                return None
                
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            return None
    
    def plot_single_result(self, filename, title=None, show=True, save_path=None):
        """绘制单个结果文件"""
        data = self.read_fortran_results(filename)
        if data is None:
            return False
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 自动生成标题
        if title is None:
            title = f"1D Convection (t={data['time']:.3f})\n"
            title += f"{data['order']}th-order {data['scheme'].upper()} + {data['rk_order']}nd-order RK"
        
        # 绘制数值解
        ax.plot(data['x'], data['numerical'], 
                label=f"Numerical ({data['scheme'].upper()}{data['order']})",
                **self.styles['numerical'])
        
        # 绘制解析解
        ax.plot(data['x'], data['analytical'],
                label="Analytical",
                **self.styles['analytical'])
        
        # 设置图形属性
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("x", fontsize=10)
        ax.set_ylabel("u", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        
        # 保存或显示
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved plot to {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()
        
        return True
        
    def format_scheme_label(self, scheme: str, order: int) -> str:
        s = scheme.upper()
        if s.startswith("WENO"):
            return f"WENO{order}"
        else:
            return f"{s}{order}"
    
    def plot_comparison(self, filenames, labels=None, title=None, show=True, save_path=None):
        """比较多个结果文件"""
        all_data = []
        print(f"plot_comparison filenames={filenames}")
        
        # 读取所有文件
        for filename in filenames:
            print(f"plot_comparison filename={filename}")
            data = self.read_fortran_results(filename)
            if data is not None:
                all_data.append(data)
        
        if not all_data:
            print("No valid data to plot")
            return False
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 子图1：所有方法比较
        for i, data in enumerate(all_data):
            #print(f"i,all_data={i,all_data}")
            if labels and i < len(labels):
                label = labels[i]
            else:
                label = f"{data['scheme'].upper()}{data['order']}"
            
            style_idx = i % len(self.styles['comparison'])
            style = self.styles['comparison'][style_idx]
            
            axes[0, 0].plot(data['x'], data['numerical'], 
                           label=label, **style)
        
        axes[0, 0].plot(all_data[0]['x'], all_data[0]['analytical'],
                       label="Analytical", **self.styles['analytical'])
        
        axes[0, 0].set_xlabel('x', fontsize=12)
        axes[0, 0].set_ylabel('u(x)', fontsize=12)
        axes[0, 0].set_title('Numerical Solutions Comparison', fontsize=14)
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # 子图2：误差分析
        for i, data in enumerate(all_data):
            if labels and i < len(labels):
                label = labels[i]
            else:
                label = f"{data['scheme'].upper()}{data['order']}"
            
            error = np.abs(data['numerical'] - data['analytical'])
            axes[0, 1].semilogy(data['x'], error, label=label)
        
        axes[0, 1].set_xlabel('x', fontsize=12)
        axes[0, 1].set_ylabel('Absolute Error', fontsize=12)
        axes[0, 1].set_title('Error Comparison', fontsize=14)
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 子图3：数值解细节（放大）
        x_min, x_max = all_data[0]['x'].min(), all_data[0]['x'].max()
        zoom_center = 1.0  # 阶跃函数位置附近
        zoom_width = 0.3
        
        for i, data in enumerate(all_data):
            if labels and i < len(labels):
                label = labels[i]
            else:
                label = f"{data['scheme'].upper()}{data['order']}"
            
            style_idx = i % len(self.styles['comparison'])
            style = self.styles['comparison'][style_idx]
            
            axes[1, 0].plot(data['x'], data['numerical'], label=label, **style)
        
        axes[1, 0].plot(all_data[0]['x'], all_data[0]['analytical'],
                       label="Analytical", **self.styles['analytical'])
        
        axes[1, 0].set_xlim(zoom_center - zoom_width/2, zoom_center + zoom_width/2)
        axes[1, 0].set_xlabel('x', fontsize=12)
        axes[1, 0].set_ylabel('u(x)', fontsize=12)
        axes[1, 0].set_title('Zoomed View (x ≈ 1.0)', fontsize=14)
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # 子图4：性能统计
        axes[1, 1].axis('off')
        
        # 计算并显示L2误差
        errors = []
        schemes = []
        
        for data in all_data:
            error = np.sqrt(np.mean((data['numerical'] - data['analytical'])**2))
            errors.append(error)
            schemes.append(f"{data['scheme'].upper()}{data['order']}")
        
        # 创建表格数据
        table_data = []
        for i, (scheme, error) in enumerate(zip(schemes, errors)):
            table_data.append([scheme, f"{error:.2e}"])
        
        # 在子图中显示表格
        table = axes[1, 1].table(cellText=table_data,
                                colLabels=['Scheme', 'L2 Error'],
                                loc='center',
                                cellLoc='center',
                                colWidths=[0.3, 0.4])
        
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        axes[1, 1].set_title('Performance Summary (L2 Error)', fontsize=14)
        
        # 设置总标题
        if title is None:
            time = all_data[0]['time']
            schemes_str = ", ".join([self.format_scheme_label(d['scheme'], d['order']) for d in all_data])
            title = f"1D Convection Comparison (t={time:.3f})\n{schemes_str}"
        
        fig.suptitle(title, fontsize=16)
        plt.tight_layout()
        
        # 保存或显示
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved comparison plot to {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()
        
        return True
        
    def get_scheme_from_filename(self, filename):
        print(f"filename={filename}")
        stem = Path(filename).stem  # e.g., "results_ENO3_40"
        parts = stem.split('_')
        if len(parts) >= 2:
            return parts[1]  # "ENO3", "WENO3", "WENO5"
        return ""        
    
    def plot_eno_weno_comparison(self, result_dir=".", save_path="eno_weno_comparison.png"):
        """自动绘制ENO/WENO对比图（类似Julia的功能）"""
        # 查找结果文件
        pattern = os.path.join(result_dir, "results_*.dat")
        files = glob.glob(pattern)
        
        if not files:
            print(f"No result files found matching {pattern}")
            return False
            
        # 重新分类
        file_info = []
        eno_files = []
        weno3_files = []
        weno5_files = []
        for f in files:
            print(f"f={f}")
            print(f"type(f)={type(f)}")
            scheme = self.get_scheme_from_filename(f)
            print(f"scheme={scheme}")
            if scheme == "ENO3":
                eno_files.append(f)
            elif scheme == "WENO3":
                weno3_files.append(f)
            elif scheme == "WENO5":
                weno5_files.append(f)            
        
        # 按求解器类型分类
        print(f"eno_files={eno_files}")
        print(f"weno3_files={weno3_files}")
        print(f"weno5_files={weno5_files}")
        
        # 选择最新的文件（如果有多组）
        files_to_plot = []
        labels = []
        
        if eno_files:
            files_to_plot.append(sorted(eno_files)[-1])
            labels.append("ENO3")
        
        if weno3_files:
            files_to_plot.append(sorted(weno3_files)[-1])
            labels.append("WENO3")
        
        if weno5_files:
            files_to_plot.append(sorted(weno5_files)[-1])
            labels.append("WENO5")
        
        if len(files_to_plot) >= 2:
            print(f"Plotting comparison of {len(files_to_plot)} solvers")
            return self.plot_comparison(files_to_plot, labels=labels, 
                                       save_path=save_path, show=True)
        else:
            print("Not enough different solvers for comparison")
            return False

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Fortran CFD Results Visualizer')
    parser.add_argument('--file', help='Plot single result file')
    parser.add_argument('--compare', nargs='+', help='Compare multiple files')
    parser.add_argument('--dir', default='.', help='Directory containing result files')
    parser.add_argument('--auto', action='store_true', help='Auto plot ENO/WENO comparison')
    parser.add_argument('--save', help='Save plot to file (without showing)')
    parser.add_argument('--no-show', action='store_true', help='Don\'t show plot')
    
    args = parser.parse_args()
    
    plotter = FortranResultsPlotter()
    
    if args.file:
        # 绘制单个文件
        plotter.plot_single_result(args.file, save_path=args.save, 
                                 show=not args.no_show)
    
    elif args.compare:
        # 比较多个文件
        plotter.plot_comparison(args.compare, save_path=args.save,
                              show=not args.no_show)
    
    elif args.auto:
        # 自动绘制比较图
        save_path = args.save or "eno_weno_comparison.png"
        plotter.plot_eno_weno_comparison(args.dir, save_path=save_path)
    
    else:
        # 默认：显示帮助
        parser.print_help()
        
        # 也显示可用的结果文件
        print("\nAvailable result files:")
        pattern = os.path.join(args.dir, "results_*.dat")
        files = glob.glob(pattern)
        
        if files:
            for f in sorted(files):
                print(f"  {os.path.basename(f)}")
            
            print(f"\nTo plot ENO/WENO comparison automatically:")
            print(f"  python plot_results.py --auto")
        else:
            print("  No result files found")

if __name__ == "__main__":
    main()