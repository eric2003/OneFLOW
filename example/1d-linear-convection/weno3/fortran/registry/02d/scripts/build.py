#!/usr/bin/env python3
"""
Fortran CFD Project Builder - 完整Python解决方案
在Python内部处理Intel oneAPI环境配置
"""

import os
import sys
import subprocess
import shutil
import argparse
import time
import platform
import tempfile
from pathlib import Path

class IntelEnvironment:
    """Intel oneAPI环境管理器"""
    
    def __init__(self):
        self.setvars_path = None
        self.env_vars = {}
        
    def find_setvars(self):
        """查找setvars.bat文件"""
        possible_paths = [
            r"C:\Program Files (x86)\Intel\oneAPI\setvars.bat",
            r"C:\Program Files\Intel\oneAPI\setvars.bat",
            r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat",
            os.path.expandvars(r"%INTEL_ONEAPI_ROOT%\setvars.bat"),
            os.path.expandvars(r"%INTEL_ONEAPI_ROOT%\compiler\latest\env\vars.bat"),
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                self.setvars_path = path
                return True
        
        return False
    
    def setup_environment(self):
        """设置Intel环境"""
        if not self.find_setvars():
            return False
        
        try:
            # 创建临时的批处理文件来捕获环境变量
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bat', delete=False) as f:
                f.write(f'@echo off\n')
                f.write(f'call "{self.setvars_path}" >nul 2>&1\n')
                f.write(f'set\n')  # 输出所有环境变量
                temp_bat = f.name
            
            # 运行批处理文件并捕获输出
            result = subprocess.run(
                [temp_bat],
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='ignore',
                shell=True
            )
            
            # 解析环境变量
            for line in result.stdout.split('\n'):
                line = line.strip()
                if '=' in line:
                    key, value = line.split('=', 1)
                    self.env_vars[key.strip()] = value.strip()
            
            # 清理临时文件
            os.unlink(temp_bat)
            
            # 更新当前进程的环境变量
            os.environ.update(self.env_vars)
            
            return True
            
        except Exception as e:
            print(f"设置Intel环境失败: {e}")
            return False
    
    def get_compiler_info(self):
        """获取编译器信息"""
        info = {}
        
        # 检查ifx编译器
        try:
            result = subprocess.run(
                ["ifx", "--version"],
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='ignore',
                env={**os.environ, **self.env_vars} if self.env_vars else os.environ
            )
            
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if 'Version' in line or '版本' in line:
                        info['ifx_version'] = line.strip()
                        break
        except:
            pass
        
        # 检查环境变量
        info['ifx_root'] = self.env_vars.get('IFX_ROOT', '')
        info['compiler_root'] = self.env_vars.get('ONEAPI_ROOT', '')
        
        return info

class BuildSystem:
    """构建系统主类"""
    
    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.build_dir = self.project_root / "build"
        self.intel_env = IntelEnvironment()
        
        # 设置控制台编码
        if sys.platform == "win32":
            try:
                import ctypes
                # 设置控制台输出为UTF-8
                ctypes.windll.kernel32.SetConsoleOutputCP(65001)
            except:
                pass
    
    def print_header(self, text):
        """打印标题"""
        print(f"\n{'='*70}")
        print(f"  {text}")
        print(f"{'='*70}\n")
    
    def print_step(self, step, total, message):
        """打印步骤"""
        print(f"[{step}/{total}] {message}...")
    
    def print_success(self, message):
        """打印成功"""
        print(f"\033[92m✓ {message}\033[0m")
    
    def print_error(self, message):
        """打印错误"""
        print(f"\033[91m✗ {message}\033[0m")
    
    def print_warning(self, message):
        """打印警告"""
        print(f"\033[93m! {message}\033[0m")
    
    def print_info(self, message):
        """打印信息"""
        print(f"\033[94mℹ {message}\033[0m")
    
    def check_prerequisites(self):
        """检查前提条件"""
        self.print_step(1, 6, "检查前提条件")
        
        # 检查Python版本
        python_version = sys.version.split()[0]
        self.print_info(f"Python版本: {python_version}")
        
        # 检查平台
        self.print_info(f"平台: {platform.system()} {platform.release()}")
        self.print_info(f"处理器核心数: {os.cpu_count()}")
        
        # 检查CMake
        try:
            result = subprocess.run(
                ["cmake", "--version"],
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='ignore'
            )
            if result.returncode == 0:
                version_line = result.stdout.split('\n')[0]
                self.print_success(f"CMake: {version_line}")
            else:
                self.print_error("CMake未找到")
                return False
        except FileNotFoundError:
            self.print_error("CMake未安装")
            return False
        
        return True
    
    def setup_intel_environment(self, args):
        """设置Intel环境"""
        self.print_step(2, 6, "配置Intel oneAPI环境")
        
        if not self.intel_env.find_setvars():
            self.print_warning("未找到Intel oneAPI setvars.bat")
            self.print_info("将尝试使用系统环境中的编译器")
            
            # 检查是否能直接访问编译器
            try:
                result = subprocess.run(
                    ["ifx", "--version"],
                    capture_output=True,
                    text=True,
                    encoding='utf-8',
                    errors='ignore'
                )
                if result.returncode == 0:
                    self.print_success("Intel编译器在系统PATH中找到")
                    return True
                else:
                    self.print_warning("Intel编译器未在PATH中找到")
            except:
                self.print_warning("无法访问Intel编译器")
            
            return True  # 继续，让CMake自己找编译器
        
        # 设置环境
        if self.intel_env.setup_environment():
            compiler_info = self.intel_env.get_compiler_info()
            
            if compiler_info.get('ifx_version'):
                self.print_success(f"Intel Fortran编译器: {compiler_info['ifx_version']}")
            elif compiler_info.get('ifx_root'):
                self.print_success(f"Intel编译器路径: {compiler_info['ifx_root']}")
            else:
                self.print_success("Intel oneAPI环境配置完成")
            
            return True
        else:
            self.print_warning("Intel环境配置失败，将继续使用系统环境")
            return True
    
    def clean_build_directory(self, args):
        """清理构建目录"""
        if args.clean and self.build_dir.exists():
            self.print_info("清理构建目录...")
            try:
                shutil.rmtree(self.build_dir)
                self.print_success("构建目录已清理")
            except Exception as e:
                self.print_error(f"清理失败: {e}")
                if not args.force:
                    return False
        return True
    
    def run_command(self, cmd, cwd=None, check=True, env=None):
        """运行命令"""
        if isinstance(cmd, list):
            cmd_str = ' '.join(str(c) for c in cmd if c)
        else:
            cmd_str = str(cmd)
        
        print(f"  \033[96m$\033[0m {cmd_str}")
        
        try:
            # 合并环境变量
            exec_env = os.environ.copy()
            if env:
                exec_env.update(env)
            if self.intel_env.env_vars:
                exec_env.update(self.intel_env.env_vars)
            
            result = subprocess.run(
                cmd,
                cwd=cwd,
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='replace',
                shell=False,
                env=exec_env
            )
            
            # 处理输出
            if result.stdout:
                for line in result.stdout.split('\n'):
                    line = line.strip()
                    if line:
                        # 高亮重要信息
                        if 'error' in line.lower():
                            print(f"    \033[91m{line}\033[0m")
                        elif 'warning' in line.lower():
                            print(f"    \033[93m{line}\033[0m")
                        elif 'success' in line.lower() or '完成' in line or '生成' in line:
                            print(f"    \033[92m{line}\033[0m")
                        else:
                            print(f"    {line}")
            
            if result.stderr:
                for line in result.stderr.split('\n'):
                    line = line.strip()
                    if line:
                        print(f"    \033[93m{line}\033[0m")
            
            if check and result.returncode != 0:
                self.print_error(f"命令执行失败，退出码: {result.returncode}")
                return False
            
            return True
            
        except Exception as e:
            self.print_error(f"命令执行异常: {e}")
            return False
    
    def configure_cmake(self, args):
        """配置CMake"""
        self.print_step(3, 6, "配置CMake项目")
        
        cmake_cmd = [
            "cmake",
            "..",
            "-G", "Visual Studio 17 2022",
            "-A", "x64",
            f"-DCMAKE_BUILD_TYPE={args.build_type}",
        ]
        
        if args.compiler == "ifx":
            cmake_cmd.extend(["-T", "fortran=ifx"])
        
        if args.verbose:
            cmake_cmd.append("-DCMAKE_VERBOSE_MAKEFILE=ON")
        
        success = self.run_command(cmake_cmd, cwd=self.build_dir, check=True)
        
        if success:
            self.print_success("CMake配置完成")
        else:
            self.print_error("CMake配置失败")
        
        return success
    
    def build_project(self, args):
        """构建项目"""
        self.print_step(4, 6, "构建项目")
        
        build_cmd = [
            "cmake",
            "--build", ".",
            "--config", args.build_type,
        ]
        
        if args.jobs > 1:
            build_cmd.append(f"-j{args.jobs}")
        
        success = self.run_command(build_cmd, cwd=self.build_dir, check=True)
        
        if success:
            self.print_success("项目构建完成")
        else:
            self.print_error("构建失败")
        
        return success
    
    def run_tests_with_environment(self, test_exe):
        """运行单个测试，确保有Intel环境"""
        try:
            # 创建临时的批处理文件来运行测试
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bat', delete=False) as f:
                if self.intel_env.setvars_path:
                    f.write(f'@echo off\n')
                    f.write(f'call "{self.intel_env.setvars_path}" >nul 2>&1\n')
                    f.write(f'"{test_exe}"\n')
                else:
                    f.write(f'@echo off\n')
                    f.write(f'"{test_exe}"\n')
                
                temp_bat = f.name
            
            # 运行测试
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
            
            return result
            
        except Exception as e:
            print(f"运行测试失败: {e}")
            return None
    
    def run_tests(self, args):
        """运行测试"""
        self.print_step(5, 6, "运行测试")
        
        # 查找测试可执行文件
        test_dir = self.build_dir / "bin" / args.build_type
        if not test_dir.exists():
            test_dir = self.build_dir / "bin"
        if not test_dir.exists():
            test_dir = self.build_dir
        
        test_files = list(test_dir.glob("test_*.exe"))
        
        if not test_files:
            self.print_warning("未找到测试程序")
            return True
        
        all_passed = True
        
        for test_exe in sorted(test_files):
            test_name = test_exe.stem
            self.print_info(f"运行测试: {test_name}")
            print(f"  {'-'*50}")
            
            # 运行测试
            result = self.run_tests_with_environment(str(test_exe))
            
            if result is None:
                self.print_error(f"  {test_name} 运行失败")
                all_passed = False
                continue
            
            # 显示输出
            if result.stdout:
                for line in result.stdout.split('\n'):
                    line = line.strip()
                    if line:
                        print(f"    {line}")
            
            if result.returncode == 0:
                self.print_success(f"  {test_name} 通过")
            else:
                self.print_error(f"  {test_name} 失败 (退出码: {result.returncode})")
                all_passed = False
                
                if result.stderr:
                    for line in result.stderr.split('\n'):
                        line = line.strip()
                        if line:
                            print(f"    \033[93m{line}\033[0m")
            
            print()  # 空行
        
        return all_passed
    
    def create_test_runner(self, args):
        """创建独立的测试运行器"""
        self.print_step(6, 6, "创建测试运行器")
        
        runner_path = self.build_dir / "run_tests.bat"
        
        content = f'''@echo off
chcp 65001 >nul

echo ========================================
echo    Fortran CFD Test Runner
echo ========================================
echo.

REM Setup Intel oneAPI environment
set "SETVARS_PATH={self.intel_env.setvars_path or ''}"
if exist "%SETVARS_PATH%" (
    call "%SETVARS_PATH%" >nul
    echo [INFO] Intel environment configured
) else (
    echo [WARNING] Intel environment not found
    echo [WARNING] Tests may fail without runtime libraries
)

echo.

REM Run all test executables
set "TEST_COUNT=0"
set "PASS_COUNT=0"

for %%f in ("bin\\{args.build_type}\\test_*.exe") do (
    set /a TEST_COUNT+=1
    echo [TEST %%f]
    echo {'-'*50}
    
    %%f
    if errorlevel 1 (
        echo [FAILED] %%f
    ) else (
        echo [PASSED] %%f
        set /a PASS_COUNT+=1
    )
    echo.
)

echo ========================================
echo Tests: %PASS_COUNT%/%TEST_COUNT% passed
if %PASS_COUNT% equ %TEST_COUNT% (
    echo [SUCCESS] All tests passed!
) else (
    echo [FAILURE] Some tests failed
)
echo ========================================

pause
'''
        
        with open(runner_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        self.print_success(f"测试运行器已创建: {runner_path}")
        self.print_info(f"使用方法: cd build && run_tests.bat")
        
        return runner_path
    
    def generate_report(self, args, build_time, tests_passed):
        """生成构建报告"""
        self.print_header("构建完成")
        
        print(f"项目: {self.project_root.name}")
        print(f"构建类型: {args.build_type}")
        print(f"编译器: {args.compiler}")
        print(f"并行作业: {args.jobs}")
        print(f"总耗时: {build_time:.1f}秒")
        print(f"测试结果: {'全部通过' if tests_passed else '有失败'}")
        
        # 显示生成的可执行文件
        bin_dir = self.build_dir / "bin" / args.build_type
        if bin_dir.exists():
            print(f"\n生成的可执行文件:")
            for exe in sorted(bin_dir.glob("*.exe")):
                size_mb = exe.stat().st_size / (1024 * 1024)
                print(f"  • {exe.name} ({size_mb:.2f} MB)")
        
        # 显示测试运行器信息
        runner_path = self.build_dir / "run_tests.bat"
        if runner_path.exists():
            print(f"\n独立测试运行器:")
            print(f"  • {runner_path.name}")
            print(f"    在Intel oneAPI环境中运行所有测试")
    
    def run(self):
        """运行构建系统"""
        parser = argparse.ArgumentParser(
            description="Fortran CFD项目构建工具",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
示例:
  %(prog)s                    # 默认构建
  %(prog)s --clean            # 清理后构建
  %(prog)s --build-type Release  # Release构建
  %(prog)s --no-tests         # 只构建，不运行测试
  %(prog)s -j8 --verbose      # 8线程并行构建，详细输出
            """
        )
        
        parser.add_argument("--build-type", choices=["Debug", "Release"], 
                          default="Debug", help="构建类型")
        parser.add_argument("--compiler", choices=["ifx", "ifort"], 
                          default="ifx", help="Fortran编译器")
        parser.add_argument("--clean", action="store_true", 
                          help="构建前清理build目录")
        parser.add_argument("--no-tests", action="store_true", 
                          help="跳过测试")
        parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count(),
                          help="并行作业数")
        parser.add_argument("--verbose", action="store_true", 
                          help="详细输出")
        parser.add_argument("--force", action="store_true",
                          help="出错时继续执行")
        
        args = parser.parse_args()
        
        # 开始构建
        start_time = time.time()
        
        self.print_header("Fortran CFD 项目构建系统")
        self.print_info(f"项目根目录: {self.project_root}")
        
        try:
            # 1. 检查前提条件
            if not self.check_prerequisites():
                if not args.force:
                    return 1
            
            # 2. 设置Intel环境
            if not self.setup_intel_environment(args):
                if not args.force:
                    return 1
            
            # 3. 清理目录
            if not self.clean_build_directory(args):
                if not args.force:
                    return 1
            
            # 确保构建目录存在
            self.build_dir.mkdir(exist_ok=True)
            
            # 4. 配置CMake
            if not self.configure_cmake(args):
                if not args.force:
                    return 1
            
            # 5. 构建项目
            if not self.build_project(args):
                if not args.force:
                    return 1
            
            # 6. 运行测试和创建测试运行器
            tests_passed = True
            if not args.no_tests:
                tests_passed = self.run_tests(args)
                
                # 创建测试运行器
                self.create_test_runner(args)  # 传递 args 参数
            
            # 7. 生成报告
            build_time = time.time() - start_time
            self.generate_report(args, build_time, tests_passed)
            
            return 0 if tests_passed else 1
            
        except KeyboardInterrupt:
            self.print_error("\n构建被用户中断")
            return 1
        except Exception as e:
            self.print_error(f"构建过程中出现未预期错误: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            return 1

def main():
    """主函数"""
    builder = BuildSystem()
    return builder.run()

if __name__ == "__main__":
    sys.exit(main())