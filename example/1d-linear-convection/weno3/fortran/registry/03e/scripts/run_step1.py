#!/usr/bin/env python3
"""
Step 1: Physics Modules Test
扩展build.py，专门用于测试物理模块
"""

import os
import sys
import subprocess
import time
from pathlib import Path

# 添加当前目录到路径，以便导入build.py的类
sys.path.insert(0, str(Path(__file__).parent))

try:
    from build import IntelEnvironment, BuildSystem
except ImportError:
    print("Error: Cannot import build.py. Make sure build.py is in the same directory.")
    sys.exit(1)

class Step1System(BuildSystem):
    """Step 1 测试系统，继承自BuildSystem"""
    
    def __init__(self):
        super().__init__()
        self.test_name = "test_physics_minimal"
        self.test_exe = None
    
    def find_test_executable(self):
        """查找测试可执行文件"""
        possible_paths = [
            self.build_dir / "bin" / "Debug" / f"{self.test_name}.exe",
            self.build_dir / "Debug" / f"{self.test_name}.exe",
            self.build_dir / f"{self.test_name}.exe",
            self.build_dir / "bin" / f"{self.test_name}.exe",
        ]
        
        for path in possible_paths:
            if path.exists():
                self.test_exe = path
                self.print_success(f"Found test executable: {path}")
                return True
        
        # 如果没有找到，尝试搜索
        self.print_warning(f"Could not find {self.test_name}.exe")
        self.print_info("Searching for test executables...")
        
        try:
            result = subprocess.run(
                ["dir", str(self.build_dir), "/s", "/b", "*.exe"],
                capture_output=True,
                text=True,
                encoding='utf-8',
                shell=True
            )
            
            if result.returncode == 0:
                test_files = [line.strip() for line in result.stdout.split('\n') 
                            if line and 'test_' in line.lower()]
                
                if test_files:
                    self.print_info("Found test files:")
                    for test_file in test_files:
                        self.print_info(f"  {test_file}")
                    return False
        except:
            pass
        
        return False
    
    def run_test_with_intel_env(self):
        """在Intel环境下运行测试"""
        if not self.test_exe:
            self.print_error("No test executable found")
            return False
        
        self.print_step(1, 2, f"Running test: {self.test_exe.name}")
        
        try:
            # 创建临时的批处理文件来运行测试（包含Intel环境）
            import tempfile
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bat', delete=False) as f:
                # 设置Intel环境
                if self.intel_env.setvars_path:
                    f.write(f'@echo off\n')
                    f.write(f'call "{self.intel_env.setvars_path}" >nul 2>&1\n')
                    f.write(f'echo [INFO] Intel environment configured\n')
                    f.write(f'echo Running test: {self.test_exe.name}\n')
                    f.write(f'echo {"-"*50}\n')
                    f.write(f'"{self.test_exe}"\n')
                else:
                    f.write(f'@echo off\n')
                    f.write(f'echo [WARNING] Intel environment not found\n')
                    f.write(f'echo Running test: {self.test_exe.name}\n')
                    f.write(f'echo {"-"*50}\n')
                    f.write(f'"{self.test_exe}"\n')
                
                temp_bat = f.name
            
            # 运行测试
            self.print_info(f"Command: {temp_bat}")
            result = subprocess.run(
                [temp_bat],
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='replace',
                shell=True
            )
            
            # 清理临时文件
            os.unlink(temp_bat)
            
            # 显示输出
            if result.stdout:
                for line in result.stdout.split('\n'):
                    line = line.strip()
                    if line:
                        # 高亮重要信息
                        if 'error' in line.lower() or 'fail' in line.lower():
                            print(f"    \033[91m{line}\033[0m")
                        elif 'warning' in line.lower():
                            print(f"    \033[93m{line}\033[0m")
                        elif 'success' in line.lower() or 'pass' in line.lower() or '✓' in line:
                            print(f"    \033[92m{line}\033[0m")
                        elif '=' in line or '---' in line:
                            print(f"    \033[96m{line}\033[0m")
                        else:
                            print(f"    {line}")
            
            if result.stderr:
                self.print_warning("Test stderr output:")
                for line in result.stderr.split('\n'):
                    line = line.strip()
                    if line:
                        print(f"    \033[93m{line}\033[0m")
            
            self.print_step(2, 2, "Test execution completed")
            
            if result.returncode == 0:
                self.print_success("Test passed")
                return True
            else:
                self.print_error(f"Test failed (exit code: {result.returncode})")
                return False
                
        except Exception as e:
            self.print_error(f"Failed to run test: {e}")
            return False
    
    def build_project_if_needed(self, args):
        """如果需要，构建项目"""
        if args.no_build:
            self.print_info("Skipping build (--no-build flag)")
            return True
        
        self.print_step(1, 3, "Building project")
        
        # 调用父类的构建方法
        build_args = argparse.Namespace()
        build_args.clean = args.clean
        build_args.build_type = "Debug"
        build_args.compiler = "ifx"
        build_args.no_tests = True  # 不运行所有测试
        build_args.jobs = os.cpu_count()
        build_args.verbose = args.verbose
        build_args.force = args.force
        
        # 清理构建目录
        if args.clean and self.build_dir.exists():
            self.print_info("Cleaning build directory...")
            import shutil
            try:
                shutil.rmtree(self.build_dir)
                self.print_success("Build directory cleaned")
            except Exception as e:
                self.print_error(f"Clean failed: {e}")
                if not args.force:
                    return False
        
        # 确保构建目录存在
        self.build_dir.mkdir(exist_ok=True)
        
        # 配置CMake
        self.print_step(2, 3, "Configuring CMake")
        cmake_cmd = [
            "cmake",
            "..",
            "-G", "Visual Studio 17 2022",
            "-A", "x64",
            "-DCMAKE_BUILD_TYPE=Debug",
            "-T", "fortran=ifx",
        ]
        
        if args.verbose:
            cmake_cmd.append("-DCMAKE_VERBOSE_MAKEFILE=ON")
        
        if not self.run_command(cmake_cmd, cwd=self.build_dir, check=not args.force):
            if not args.force:
                return False
        
        # 构建项目
        self.print_step(3, 3, "Building project")
        build_cmd = [
            "cmake",
            "--build", ".",
            "--config", "Debug",
        ]
        
        if args.jobs > 1:
            build_cmd.append(f"-j{args.jobs}")
        
        if not self.run_command(build_cmd, cwd=self.build_dir, check=not args.force):
            if not args.force:
                return False
        
        self.print_success("Build completed")
        return True
    
    def run(self):
        """运行Step 1测试"""
        parser = argparse.ArgumentParser(
            description="Step 1: Physics Modules Test",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
示例:
  %(prog)s                    # 默认运行
  %(prog)s --clean            # 清理后构建并测试
  %(prog)s --no-build         # 只运行测试，不重新构建
  %(prog)s --verbose          # 详细输出
  %(prog)s -j4               # 使用4个并行作业构建
            """
        )
        
        parser.add_argument("--clean", action="store_true", 
                          help="构建前清理build目录")
        parser.add_argument("--no-build", action="store_true",
                          help="不重新构建，直接运行测试")
        parser.add_argument("--verbose", action="store_true", 
                          help="详细输出")
        parser.add_argument("--force", action="store_true",
                          help="出错时继续执行")
        parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count(),
                          help="并行作业数")
        
        args = parser.parse_args()
        
        # 开始测试
        start_time = time.time()
        
        self.print_header("Step 1: Physics Modules Implementation Test")
        self.print_info(f"项目根目录: {self.project_root}")
        
        try:
            # 1. 设置Intel环境
            self.print_step(1, 4, "Setting up Intel oneAPI environment")
            if not self.setup_intel_environment(args):
                if not args.force:
                    self.print_error("Intel environment setup failed")
                    return 1
                self.print_warning("Intel environment setup failed, continuing...")
            
            # 2. 构建项目（如果需要）
            self.print_step(2, 4, "Building project if needed")
            if not self.build_project_if_needed(args):
                if not args.force:
                    return 1
            
            # 3. 查找测试可执行文件
            self.print_step(3, 4, "Finding test executable")
            if not self.find_test_executable():
                if not args.force:
                    return 1
                self.print_warning("Test executable not found, but continuing due to --force")
                return 0
            
            # 4. 运行测试
            self.print_step(4, 4, "Running physics module test")
            test_passed = self.run_test_with_intel_env()
            
            # 生成报告
            test_time = time.time() - start_time
            self.print_header("Step 1 Complete")
            
            print(f"测试: {'通过 ✓' if test_passed else '失败 ✗'}")
            print(f"测试程序: {self.test_exe.name if self.test_exe else '未找到'}")
            print(f"总耗时: {test_time:.1f}秒")
            
            if test_passed:
                print(f"\n下一步: 更新配置以包含物理设置")
                print(f"建议: 修改config.f90，添加physics相关字段")
                return 0
            else:
                if args.force:
                    self.print_warning("测试失败，但由于--force标志继续执行")
                    return 0
                return 1
            
        except KeyboardInterrupt:
            self.print_error("\n测试被用户中断")
            return 1
        except Exception as e:
            self.print_error(f"测试过程中出现未预期错误: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            return 1

def main():
    """主函数"""
    system = Step1System()
    return system.run()

if __name__ == "__main__":
    # 需要导入argparse
    import argparse
    sys.exit(main())