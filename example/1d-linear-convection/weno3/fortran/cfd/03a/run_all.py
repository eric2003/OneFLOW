#!/usr/bin/env python3
"""
OneFLOW-CFD Complete Automation System
Author: Your Name
Date: 2024
Description: Automated build, run, and visualization pipeline
"""

import os
import sys
import subprocess
import shutil
import time
from pathlib import Path
import argparse

class CFDAutomation:
    """自动化CFD求解器的完整流程"""
    
    def __init__(self, build_type="Release", verbose=False):
        self.project_root = Path.cwd()
        self.build_dir = self.project_root / "build"
        self.bin_dir = None
        self.executable = None
        self.build_type = build_type
        self.verbose = verbose
        self.results_dir = self.project_root / "results"
        # 初始化combined_env，避免属性不存在报错
        self.combined_env = os.environ.copy()
        
    def print_banner(self):
        """打印漂亮的横幅"""
        print("\n" + "="*70)
        print("    OneFLOW-CFD Automation System")
        print("    1D Convection: ENO3 vs WENO3 Comparison")
        print("="*70 + "\n")
        
    def print_step(self, step_num, description):
        """打印步骤信息"""
        print(f"[{step_num}/7] {description}")
        print("-" * 50)
        
    def check_environment(self):
        """检查运行环境"""
        self.print_step(2, "Checking Environment")  # 调整步骤号，匹配后续流程
        
        # 检查Python
        try:
            import numpy, matplotlib
            print("✓ NumPy and Matplotlib are available")
        except ImportError as e:
            print(f"✗ Missing Python packages: {e}")
            print("Installing required packages...")
            try:
                subprocess.run([sys.executable, "-m", "pip", "install", "numpy", "matplotlib", "-q"],
                              capture_output=True,
                              text=True,
                              encoding='utf-8',
                              errors='ignore',
                              shell=True)
                print("✓ Packages installed")
            except Exception as install_e:
                print(f"✗ Package installation failed: {install_e}")
                sys.exit(1)
                
        # 检查Fortran编译器
        compilers = ["ifx", "ifort", "gfortran", "flang"]
        found = False
        for compiler in compilers:
            try:
                result = subprocess.run([compiler, "--version"], 
                                      env=self.combined_env,
                                      capture_output=True,
                                      text=True,
                                      encoding='utf-8',  # 替换为utf-8，避免gbk解码失败
                                      errors='ignore',    # 忽略无效字符，兜底报错
                                      shell=True)
                if result.returncode == 0:
                    print(f"✓ Found Fortran compiler: {compiler}")
                    found = True
                    break
            except (FileNotFoundError, Exception):
                continue
                
        if not found:
            print("⚠ Warning: No Fortran compiler found in PATH")
            print("  On Windows, please run from Intel oneAPI command prompt")
            print("  On Linux/Mac, install: ifx, gfortran, or flang")
            
        # 检查CMake
        try:
            subprocess.run(["cmake", "--version"], 
                           env=self.combined_env,
                           capture_output=True,
                           text=True,
                           encoding='utf-8',  # 替换为utf-8
                           errors='ignore',    # 兜底无效字符
                           shell=True)            
            print("✓ CMake is available")
        except FileNotFoundError:
            print("✗ CMake not found. Please install CMake 3.10+")
            sys.exit(1)
            
    def setup_intel_environment(self):
        """一次性提取setvars.bat配置后的所有环境变量（仅执行一次）"""
        self.print_step(1, "Setting up Intel oneAPI Environment")  # 明确步骤1
        
        cmd_command = r'cmd.exe /C ""C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && set"'

        try:
            result = subprocess.run(
                cmd_command,
                capture_output=True,
                text=True,
                encoding='gbk',  # 保留gbk解析cmd输出
                errors='ignore', # 新增：忽略无法解码的字符，避免报错
                shell=True,
                timeout=30
            )
        except Exception as e:
            raise Exception(f"提取oneAPI环境失败：{str(e)}")

        if result.returncode != 0:
            raise Exception(f"setvars.bat执行失败：{result.stderr}")

        # 解析环境变量（忽略非有效行）
        env_dict = {}
        for line in result.stdout.splitlines():
            line = line.strip()
            if '=' in line and not line.startswith('='):
                try:
                    key, value = line.split('=', 1)
                    env_dict[key.strip()] = value.strip()
                except:
                    continue  # 忽略分割失败的行，避免中断解析
                
        oneapi_env = env_dict

        # 合并环境变量（系统原有环境 + oneAPI环境）
        self.combined_env = os.environ.copy()
        self.combined_env.update(oneapi_env)                

        print("✓ Intel oneAPI environment loaded successfully")
        return oneapi_env
            
    def configure_cmake(self):
        """配置CMake"""
        self.print_step(3, "Configuring with CMake")
        
        # 清理旧的构建目录
        if self.build_dir.exists():
            print("Cleaning previous build...")
            shutil.rmtree(self.build_dir)
            
        # 创建构建目录
        self.build_dir.mkdir(exist_ok=True)
        
        # CMake命令
        cmake_cmd = ["cmake", "-S", str(self.project_root), "-B", str(self.build_dir)]
        
        # 平台特定选项
        if sys.platform == "win32":
            cmake_cmd.extend(["-G", "Visual Studio 17 2022", "-A", "x64"])
            cmake_cmd.extend(["-T", "fortran=ifx"])
        else:
            cmake_cmd.extend(["-DCMAKE_BUILD_TYPE=" + self.build_type])
            
        print(f"Running: {' '.join(cmake_cmd)}")
        
        try:
            # 核心修改：指定编码为utf-8，添加errors='ignore'
            result = subprocess.run(cmake_cmd, check=True, 
                                  capture_output=not self.verbose,
                                  env=self.combined_env,
                                  text=True,
                                  encoding='utf-8',
                                  errors='ignore',
                                  shell=True) 
            if self.verbose and result.stdout:
                print(result.stdout)
            print("✓ CMake configuration successful")
        except subprocess.CalledProcessError as e:
            print(f"✗ CMake configuration failed:")
            if e.stderr:
                print(e.stderr)
            sys.exit(1)
            
    def build_project(self):
        """构建项目"""
        self.print_step(4, "Building Project")
        
        build_cmd = ["cmake", "--build", str(self.build_dir), "--config", self.build_type]
        
        if self.verbose:
            build_cmd.append("--verbose")
            
        print(f"Running: {' '.join(build_cmd)}")
        
        try:
            # 核心修改：指定编码为utf-8，添加errors='ignore'，解决UnicodeDecodeError
            result = subprocess.run(build_cmd, check=True,
                                  capture_output=not self.verbose,
                                  env=self.combined_env,
                                  text=True,
                                  encoding='utf-8',
                                  errors='ignore',
                                  shell=True)   
            if self.verbose and result.stdout:
                print(result.stdout)
            print("✓ Build successful")
            
            # 查找可执行文件
            self._find_executable()
            
        except subprocess.CalledProcessError as e:
            print(f"✗ Build failed:")
            if e.stderr:
                print(e.stderr)
            sys.exit(1)
            
    def _find_executable(self):
        """查找生成的可执行文件"""
        search_paths = []
        
        if sys.platform == "win32":
            search_paths = [
                self.build_dir / "bin" / self.build_type / "oneflow_cfd.exe",
                self.build_dir / self.build_type / "oneflow_cfd.exe",
            ]
        else:
            search_paths = [
                self.build_dir / "oneflow_cfd",
                self.build_dir / "bin" / "oneflow_cfd",
            ]
            
        for path in search_paths:
            if path.exists():
                self.executable = path
                self.bin_dir = path.parent
                print(f"✓ Found executable: {self.executable}")
                return
                
        print("⚠ Could not find executable automatically")
        print("  Searched in:", [str(p) for p in search_paths])
        
    def run_simulation(self):
        """运行CFD模拟"""
        self.print_step(5, "Running CFD Simulation")
        
        if not self.executable or not self.executable.exists():
            print("✗ Executable not found!")
            return False
            
        # 确保postprocess.py在运行目录
        postprocess_src = self.project_root / "postprocess.py"
        if postprocess_src.exists():
            shutil.copy2(postprocess_src, self.bin_dir / "postprocess.py")
            
        # 运行可执行文件
        print(f"Running: {self.executable.name}")
        print("-" * 50)
        
        try:
            original_dir = os.getcwd()
            os.chdir(self.bin_dir)
            
            result = subprocess.run([str(self.executable)], 
                                  check=True,
                                  capture_output=True,
                                  text=True,
                                  encoding='utf-8',
                                  errors='ignore',
                                  env=self.combined_env,
                                  shell=True)   
            
            # 打印输出
            if result.stdout:
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    if any(keyword in line.lower() for keyword in 
                          ['error', 'warning', 'failed', 'success']):
                        if 'error' in line.lower() or 'failed' in line.lower():
                            print(f"✗ {line}")
                        elif 'warning' in line.lower():
                            print(f"⚠ {line}")
                        else:
                            print(f"✓ {line}")
                    else:
                        print(f"  {line}")
                        
            if result.stderr:
                print("\nStderr output:")
                print(result.stderr)
                
            print("-" * 50)
            print("✓ Simulation completed")
            
            os.chdir(original_dir)
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"✗ Simulation failed with exit code {e.returncode}")
            if e.stdout:
                print("Output:", e.stdout)
            os.chdir(original_dir)
            return False
            
    def generate_visualization(self):
        """生成可视化结果"""
        self.print_step(6, "Generating Visualization")
        
        original_dir = os.getcwd()
        os.chdir(self.bin_dir)        
        
        postprocess_script = self.bin_dir / "postprocess.py"
        if not postprocess_script.exists():
            print("✗ postprocess.py not found in output directory")
            os.chdir(original_dir)
            return False
            
        print(f"Running: python {postprocess_script.name}")
        
        try:
            result = subprocess.run([sys.executable, str(postprocess_script)],
                                  check=True,
                                  capture_output=True,
                                  text=True,
                                  encoding='utf-8',  # 新增：指定utf-8编码
                                  errors='ignore',    # 新增：忽略无效字符
                                  env=self.combined_env,
                                  shell=True)   

            if result.stdout and self.verbose:
                print(result.stdout)
                
            print("✓ Visualization generated")
            os.chdir(original_dir)
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"✗ Visualization failed: {e}")
            os.chdir(original_dir)
            return False
            
    def organize_results(self):
        """整理结果文件"""
        self.print_step(7, "Organizing Results")
        
        # 创建结果目录
        self.results_dir.mkdir(exist_ok=True)
        
        # 复制结果文件
        result_patterns = ["*results.txt", "comparison.png", "*.dat"]
        
        for pattern in result_patterns:
            for result_file in self.bin_dir.glob(pattern):
                if result_file.is_file():
                    dest = self.results_dir / result_file.name
                    # 避免同名文件覆盖提示，直接复制
                    shutil.copy2(result_file, dest)
                    print(f"✓ Copied: {result_file.name}")
                    
        # 创建结果摘要
        self._create_summary()
        
        print(f"\n✓ Results organized in: {self.results_dir}")
        
    def _create_summary(self):
        """创建结果摘要"""
        summary_file = self.results_dir / "simulation_summary.txt"
        
        with open(summary_file, 'w', encoding='utf-8') as f:  # 指定utf-8写入，避免中文乱码
            f.write("OneFLOW-CFD Simulation Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Build Type: {self.build_type}\n")
            f.write(f"Platform: {sys.platform}\n")
            f.write("\nFiles Generated:\n")
            
            for file in sorted(self.results_dir.glob("*")):
                if file.is_file():
                    size_kb = file.stat().st_size / 1024
                    f.write(f"  - {file.name:20} {size_kb:6.1f} KB\n")
                    
    def run_complete_pipeline(self):
        """运行完整管道"""
        self.print_banner()
        
        steps = [
            self.setup_intel_environment,
            self.check_environment,            
            self.configure_cmake,
            self.build_project,
            self.run_simulation,
            self.generate_visualization,
            self.organize_results,
        ]
        
        for i, step in enumerate(steps, 1):
            try:
                step()
            except Exception as e:
                print(f"\n✗ Step {i} failed: {e}")
                print("\nDebug information:")
                print(f"  Project root: {self.project_root}")
                print(f"  Build dir: {self.build_dir}")
                print(f"  Executable: {self.executable}")
                sys.exit(1)
                
        self._print_completion_message()
        
    def _print_completion_message(self):
        """打印完成信息"""
        print("\n" + "="*70)
        print("    OneFLOW-CFD Automation Complete!")
        print("="*70)
        print("\nGenerated Files:")
        print("  Results Directory: results/")
        
        for file in sorted(self.results_dir.glob("*")):
            if file.is_file():
                print(f"    - {file.name}")
                
        print("\nTo view the plot:")
        plot_file = self.results_dir / "comparison.png"
        if plot_file.exists():
            if sys.platform == "win32":
                print(f"    start results/comparison.png")
            elif sys.platform == "darwin":
                print(f"    open results/comparison.png")
            else:
                print(f"    xdg-open results/comparison.png")
                
        print("\n" + "="*70)

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="OneFLOW-CFD Automation System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Run complete pipeline
  %(prog)s --verbose          # Verbose output
  %(prog)s --build Debug      # Debug build
  %(prog)s --clean            # Clean before building
        """
    )
    
    parser.add_argument("--build", choices=["Release", "Debug", "RelWithDebInfo"],
                       default="Release", help="Build type")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    parser.add_argument("--clean", "-c", action="store_true",
                       help="Clean build directory before building")
    parser.add_argument("--no-vis", action="store_true",
                       help="Skip visualization")
    parser.add_argument("--no-run", action="store_true",
                       help="Skip running simulation")
    
    args = parser.parse_args()
    
    # 创建自动化实例
    automator = CFDAutomation(build_type=args.build, verbose=args.verbose)
    
    # 清理选项
    if args.clean and automator.build_dir.exists():
        print("Cleaning build directory...")
        shutil.rmtree(automator.build_dir)
        
    # 运行管道
    try:
        automator.run_complete_pipeline()
    except KeyboardInterrupt:
        print("\n\n⚠ Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\n✗ Automation failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()