#!/usr/bin/env python3
"""
Fortran CFD Project Builder (Python)
支持 Intel oneAPI 环境的自动化构建脚本
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple, Optional
import platform

class Color:
    """终端颜色代码"""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_step(step_num: int, total_steps: int, message: str):
    """打印步骤信息"""
    print(f"{Color.CYAN}[{step_num}/{total_steps}]{Color.END} {message}...")

def print_success(message: str):
    """打印成功信息"""
    print(f"{Color.GREEN}✓{Color.END} {message}")

def print_error(message: str):
    """打印错误信息"""
    print(f"{Color.RED}[ERROR]{Color.END} {message}")

def print_warning(message: str):
    """打印警告信息"""
    print(f"{Color.YELLOW}[WARN]{Color.END} {message}")

def run_command(cmd: List[str], cwd: Optional[str] = None, check: bool = True) -> Tuple[int, str]:
    """运行命令并返回结果"""
    print(f"{Color.GRAY}Running: {' '.join(cmd)}{Color.END}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )
        
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(f"{Color.YELLOW}{result.stderr}{Color.END}")
            
        if check and result.returncode != 0:
            print_error(f"Command failed with exit code: {result.returncode}")
            
        return result.returncode, result.stdout
        
    except FileNotFoundError as e:
        print_error(f"Command not found: {cmd[0]}")
        if check:
            raise
        return 1, str(e)

def setup_intel_environment():
    """设置 Intel oneAPI 环境"""
    setvars_path = r"C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
    
    if not os.path.exists(setvars_path):
        print_error(f"Intel oneAPI setvars.bat not found at: {setvars_path}")
        print("Please check your Intel oneAPI installation.")
        return False
    
    # 创建临时的 batch 文件来设置环境
    temp_bat = "setup_intel_env.bat"
    with open(temp_bat, 'w') as f:
        f.write(f'@echo off\n')
        f.write(f'call "{setvars_path}" > nul 2>&1\n')
        f.write(f'set\n')  # 输出所有环境变量
    
    try:
        # 运行 batch 文件并捕获环境变量
        env_result = subprocess.run(
            [temp_bat],
            capture_output=True,
            text=True,
            shell=True,
            encoding='utf-8',
            errors='ignore'
        )
        
        # 解析环境变量并设置
        for line in env_result.stdout.split('\n'):
            if '=' in line:
                key, value = line.split('=', 1)
                os.environ[key.strip()] = value.strip()
        
        os.remove(temp_bat)
        print_success("Intel oneAPI environment set")
        return True
        
    except Exception as e:
        print_error(f"Failed to set Intel environment: {e}")
        if os.path.exists(temp_bat):
            os.remove(temp_bat)
        return False

def main():
    """主函数"""
    print(f"{Color.CYAN}{'='*50}{Color.END}")
    print(f"{Color.BOLD}  Fortran CFD Project Builder (Python){Color.END}")
    print(f"{Color.CYAN}{'='*50}{Color.END}\n")
    
    total_steps = 5
    
    # 步骤1: 设置 Intel 环境
    print_step(1, total_steps, "Setting up Intel oneAPI environment")
    if not setup_intel_environment():
        input("Press Enter to exit...")
        sys.exit(1)
    
    # 步骤2: 清理构建目录
    print_step(2, total_steps, "Cleaning build directory")
    build_dir = Path("build")
    if build_dir.exists():
        try:
            shutil.rmtree(build_dir)
            print_success("Old build directory removed")
        except Exception as e:
            print_error(f"Failed to remove build directory: {e}")
            sys.exit(1)
    else:
        print("- No existing build directory")
    
    build_dir.mkdir(exist_ok=True)
    os.chdir(build_dir)
    
    # 步骤3: 配置项目
    print_step(3, total_steps, "Configuring project with Intel Fortran")
    cmake_cmd = [
        "cmake",
        "-G", "Visual Studio 17 2022",
        "-A", "x64",
        "-T", "fortran=ifx",
        ".."
    ]
    
    return_code, _ = run_command(cmake_cmd, check=False)
    if return_code != 0:
        print_error("CMake configuration failed")
        input("Press Enter to exit...")
        sys.exit(1)
    print_success("CMake configuration successful")
    
    # 步骤4: 构建项目
    print_step(4, total_steps, "Building project")
    build_cmd = ["cmake", "--build", ".", "--config", "Debug"]
    
    return_code, _ = run_command(build_cmd, check=False)
    if return_code != 0:
        print_error("Build failed")
        input("Press Enter to exit...")
        sys.exit(1)
    print_success("Build successful")
    
    # 步骤5: 运行测试
    print_step(5, total_steps, "Running tests")
    print(f"{Color.CYAN}{'='*50}{Color.END}")
    
    test_results = []
    
    # 运行简单测试
    test_simple = Path("Debug/test_simple.exe")
    if test_simple.exists():
        print(f"{Color.MAGENTA}[TEST] Simple functionality test...{Color.END}")
        print(f"{Color.DARK_GRAY}{'-'*50}{Color.END}")
        return_code, output = run_command([str(test_simple)], check=False)
        if return_code == 0:
            print_success("Simple test passed")
            test_results.append(("Simple Test", "PASSED", Color.GREEN))
        else:
            print_error("Simple test failed")
            test_results.append(("Simple Test", "FAILED", Color.RED))
        print()
    else:
        print_warning("test_simple.exe not found")
    
    # 运行工厂测试
    test_factory = Path("Debug/test_factory.exe")
    if test_factory.exists():
        print(f"{Color.MAGENTA}[TEST] Factory pattern test...{Color.END}")
        print(f"{Color.DARK_GRAY}{'-'*50}{Color.END}")
        return_code, output = run_command([str(test_factory)], check=False)
        if return_code == 0:
            print_success("Factory test passed")
            test_results.append(("Factory Test", "PASSED", Color.GREEN))
        else:
            print_error("Factory test failed")
            test_results.append(("Factory Test", "FAILED", Color.RED))
        print()
    else:
        print_warning("test_factory.exe not found")
    
    # 显示测试总结
    print(f"{Color.CYAN}{'='*50}{Color.END}")
    print(f"{Color.BOLD}TEST SUMMARY:{Color.END}")
    print(f"{Color.DARK_GRAY}{'-'*50}{Color.END}")
    
    for test_name, status, color in test_results:
        print(f"{color}{test_name}: {status}{Color.END}")
    
    print(f"{Color.CYAN}{'='*50}{Color.END}")
    print(f"{Color.BOLD}Build directory:{Color.END} {os.getcwd()}")
    print(f"{Color.CYAN}{'='*50}{Color.END}")
    
    # 如果有测试失败，等待用户确认
    if any(status == "FAILED" for _, status, _ in test_results):
        input("Press Enter to exit...")

if __name__ == "__main__":
    # 添加颜色支持
    if platform.system() == "Windows":
        # Windows 需要启用 ANSI 转义序列
        os.system("")
    
    main()