#!/usr/bin/env python3
"""
Fortran CFD 项目构建脚本 (Python版)
支持 Intel oneAPI 环境的自动化构建
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path
import argparse
import platform
import time

class Colors:
    """终端颜色"""
    if platform.system() == "Windows":
        # Windows 启用 ANSI
        os.system("")
    
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    MAGENTA = '\033[95m'
    DARK_GRAY = '\033[90m'  # 添加这个
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'    

def print_color(text, color=Colors.END):
    """彩色打印"""
    print(f"{color}{text}{Colors.END}")

def print_header(text):
    """打印标题"""
    print_color("\n" + "="*60, Colors.CYAN)
    print_color(f"  {text}", Colors.BOLD + Colors.CYAN)
    print_color("="*60 + "\n", Colors.CYAN)

def print_step(step, total, message):
    """打印步骤"""
    print_color(f"[{step}/{total}] {message}...", Colors.YELLOW)

def print_success(message):
    """打印成功"""
    print_color(f"✓ {message}", Colors.GREEN)

def print_error(message):
    """打印错误"""
    print_color(f"✗ {message}", Colors.RED)

def print_warning(message):
    """打印警告"""
    print_color(f"! {message}", Colors.YELLOW)

def run_command(cmd, cwd=None, check=True, capture=True):
    """运行命令"""
    print_color(f"  $ {' '.join(cmd)}", Colors.BLUE)
    
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=capture,
            text=True,
            encoding='utf-8',
            errors='ignore',
            shell=False
        )
        
        if capture and result.stdout:
            print(result.stdout)
        if capture and result.stderr:
            print_color(result.stderr, Colors.YELLOW)
        
        if check and result.returncode != 0:
            print_error(f"Command failed with exit code: {result.returncode}")
            return False
        
        return True
        
    except FileNotFoundError as e:
        print_error(f"Command not found: {cmd[0]}")
        if check:
            raise
        return False
    except Exception as e:
        print_error(f"Command execution failed: {e}")
        return False

def setup_intel_environment():
    """设置 Intel oneAPI 环境"""
    setvars_paths = [
        r"C:\Program Files (x86)\Intel\oneAPI\setvars.bat",
        r"C:\Program Files\Intel\oneAPI\setvars.bat",
    ]
    
    setvars_path = None
    for path in setvars_paths:
        if os.path.exists(path):
            setvars_path = path
            break
    
    if not setvars_path:
        print_error("Intel oneAPI setvars.bat not found.")
        print_warning("Please install Intel oneAPI or update the path in build.py")
        return False
    
    # 创建临时的 batch 文件
    temp_bat = "temp_setvars.bat"
    with open(temp_bat, 'w') as f:
        f.write(f'@echo off\n')
        f.write(f'call "{setvars_path}" > nul 2>&1\n')
        f.write(f'set\n')
    
    try:
        # 运行并捕获环境变量
        result = subprocess.run(
            [temp_bat],
            capture_output=True,
            text=True,
            shell=True,
            encoding='utf-8',
            errors='ignore'
        )
        
        # 解析并设置环境变量
        for line in result.stdout.split('\n'):
            if '=' in line:
                key, value = line.split('=', 1)
                os.environ[key.strip()] = value.strip()
        
        os.remove(temp_bat)
        print_success("Intel oneAPI environment configured")
        return True
        
    except Exception as e:
        print_error(f"Failed to setup Intel environment: {e}")
        if os.path.exists(temp_bat):
            os.remove(temp_bat)
        return False

def build_project(args):
    """构建项目主函数"""
    start_time = time.time()
    
    print_header(f"Fortran CFD Project Builder")
    print_color(f"Build type: {args.build_type}", Colors.CYAN)
    print_color(f"Run tests: {args.run_tests}", Colors.CYAN)
    print()
    
    # 获取项目根目录（脚本所在目录的父目录）
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    os.chdir(project_root)
    
    print_color(f"Project root: {project_root}", Colors.BLUE)
    
    # 步骤1: 设置 Intel 环境
    print_step(1, 4, "Setting up Intel Fortran compiler")
    if not setup_intel_environment():
        return False
    
    # 步骤2: 准备构建目录
    print_step(2, 4, "Preparing build directory")
    build_dir = project_root / "build"
    
    if args.clean and build_dir.exists():
        try:
            shutil.rmtree(build_dir)
            print_success("Cleaned build directory")
        except Exception as e:
            print_error(f"Failed to clean build directory: {e}")
            if not args.force:
                return False
    
    build_dir.mkdir(exist_ok=True)
    os.chdir(build_dir)
    
    # 步骤3: 配置项目
    print_step(3, 4, "Configuring project")
    
    cmake_cmd = [
        "cmake",
        "-G", "Visual Studio 17 2022",
        "-A", "x64",
        "-T", "fortran=ifx",
        f"-DCMAKE_BUILD_TYPE={args.build_type}",
        ".."
    ]
    
    if not run_command(cmake_cmd, check=not args.force):
        if not args.force:
            return False
    
    # 步骤4: 构建项目
    print_step(4, 4, "Building project")
    
    build_cmd = [
        "cmake",
        "--build", ".",
        "--config", args.build_type,
        f"-j{args.jobs}" if args.jobs > 1 else ""
    ]
    build_cmd = [c for c in build_cmd if c]  # 移除空字符串
    
    if not run_command(build_cmd, check=not args.force):
        if not args.force:
            return False
    
    build_time = time.time() - start_time
    print_success(f"Build completed in {build_time:.1f} seconds")
    
    # 运行测试
    if args.run_tests:
        print_header("Running Tests")
        run_tests(args.build_type)
    
    print_header("Build Summary")
    print_color(f"Build directory: {build_dir}", Colors.GREEN)
    print_color(f"Build type: {args.build_type}", Colors.GREEN)
    print_color(f"Total time: {build_time:.1f}s", Colors.GREEN)
    
    return True

def run_tests(build_type):
    """运行测试"""
    tests = [
        ("test_simple", "Simple functionality test"),
        ("test_factory", "Factory pattern test"),
    ]
    
    for test_name, description in tests:
        # 修复：在这里直接访问当前循环的 test_name
        _run_single_test(test_name, description, build_type)

def _run_single_test(test_name, description, build_type):
    """运行单个测试（修复作用域问题）"""
    # 可能的测试可执行文件路径
    possible_paths = [
        Path(f"./{build_type}/{test_name}.exe"),
        Path(f"./tests/{build_type}/{test_name}.exe"),
        Path(f"./tests/{test_name}.exe"),
        Path(f"{test_name}.exe"),
        Path(f"./{test_name}.exe"),
    ]
    
    test_exe = None
    
    # 尝试所有可能的路径
    for path in possible_paths:
        if path.exists():
            test_exe = path
            break
    
    if test_exe:
        print_color(f"\n[TEST] {description}...", Colors.MAGENTA)
        print_color("-" * 50, Colors.DARK_GRAY)
        
        try:
            result = subprocess.run(
                [str(test_exe)],
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='ignore'
            )
            
            if result.stdout:
                print(result.stdout)
            
            if result.returncode == 0:
                print_success(f"{test_name} passed")
            else:
                print_error(f"{test_name} failed (exit code: {result.returncode})")
                if result.stderr:
                    print_color(result.stderr, Colors.YELLOW)
                    
        except Exception as e:
            print_error(f"{test_name} failed: {e}")
    else:
        print_warning(f"{test_name}.exe not found. Searched paths:")
        for path in possible_paths:
            print_color(f"  - {path}", Colors.YELLOW)

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Fortran CFD Project Builder")
    
    parser.add_argument(
        "--build-type",
        choices=["Debug", "Release", "RelWithDebInfo", "MinSizeRel"],
        default="Debug",
        help="Build type (default: Debug)"
    )
    
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Clean build directory before building"
    )
    
    parser.add_argument(
        "--no-tests",
        action="store_true",
        help="Skip running tests"
    )
    
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=0,
        help="Number of parallel jobs (0 = use all cores)"
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Continue on errors"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    if args.jobs == 0:
        import multiprocessing
        args.jobs = multiprocessing.cpu_count()
    
    args.run_tests = not args.no_tests
    
    try:
        success = build_project(args)
        if not success and not args.force:
            sys.exit(1)
    except KeyboardInterrupt:
        print_color("\nBuild interrupted by user", Colors.YELLOW)
        sys.exit(1)
    except Exception as e:
        print_error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()