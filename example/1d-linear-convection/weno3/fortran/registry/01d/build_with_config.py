#!/usr/bin/env python3
"""
配置文件驱动的构建脚本
"""

import yaml
import os
import sys
from pathlib import Path
from build import (  # 复用上面的 build.py 中的函数
    Color, print_step, print_success, 
    print_error, print_warning, run_command
)

def load_config(config_file="build_config.yaml"):
    """加载配置文件"""
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except FileNotFoundError:
        print_error(f"Config file not found: {config_file}")
        return None
    except yaml.YAMLError as e:
        print_error(f"Error parsing config file: {e}")
        return None

def setup_intel_from_config(config):
    """根据配置设置 Intel 环境"""
    setvars_path = config['intel']['setvars_path']
    
    if not os.path.exists(setvars_path):
        print_error(f"Intel setvars.bat not found: {setvars_path}")
        return False
    
    # 这里可以复用 build.py 中的 setup_intel_environment 函数
    # 或者简化的版本
    os.environ['PATH'] = f"C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\2025.2\\bin\\intel64;{os.environ['PATH']}"
    return True

def main():
    """主函数 - 配置文件版本"""
    config = load_config()
    if not config:
        sys.exit(1)
    
    print(f"{Color.CYAN}{'='*60}{Color.END}")
    print(f"{Color.BOLD}  {config['project']['name']} v{config['project']['version']} - Builder{Color.END}")
    print(f"{Color.CYAN}{'='*60}{Color.END}\n")
    
    total_steps = 4
    
    # 步骤1: 设置环境
    print_step(1, total_steps, f"Setting up {config['intel']['compiler']} compiler")
    if not setup_intel_from_config(config):
        input("Press Enter to exit...")
        sys.exit(1)
    print_success(f"{config['intel']['compiler']} environment set")
    
    # 步骤2: 准备构建目录
    print_step(2, total_steps, "Preparing build directory")
    build_dir = Path("build")
    
    if config['build']['clean_before_build'] and build_dir.exists():
        import shutil
        try:
            shutil.rmtree(build_dir)
            print_success("Cleaned build directory")
        except Exception as e:
            print_error(f"Failed to clean: {e}")
            sys.exit(1)
    
    build_dir.mkdir(exist_ok=True)
    os.chdir(build_dir)
    
    # 步骤3: 配置和构建
    print_step(3, total_steps, "Configuring and building")
    
    cmake_cmd = [
        "cmake",
        "-G", config['cmake']['generator'],
        "-A", config['cmake']['platform'],
        "-T", config['cmake']['toolset'],
        config['cmake']['source_dir']
    ]
    
    return_code, _ = run_command(cmake_cmd, check=False)
    if return_code != 0:
        print_error("CMake configuration failed")
        sys.exit(1)
    
    build_cmd = [
        "cmake", "--build", ".", 
        "--config", config['cmake']['build_type']
    ]
    
    return_code, _ = run_command(build_cmd, check=False)
    if return_code != 0:
        print_error("Build failed")
        sys.exit(1)
    
    print_success("Build completed")
    
    # 步骤4: 运行测试
    if config['build']['run_tests']:
        print_step(4, total_steps, "Running tests")
        print(f"{Color.CYAN}{'-'*60}{Color.END}")
        
        for test in config['tests']:
            test_exe = Path(test['executable'])
            if test_exe.exists():
                print(f"{Color.MAGENTA}[TEST] {test['description']}...{Color.END}")
                return_code, output = run_command([str(test_exe)], check=False)
                if return_code == 0:
                    print_success(f"{test['name']} passed")
                else:
                    print_error(f"{test['name']} failed")
                print()
            else:
                print_warning(f"{test['name']} executable not found: {test_exe}")
    
    print(f"{Color.CYAN}{'='*60}{Color.END}")
    print(f"{Color.BOLD}Build completed successfully!{Color.END}")
    print(f"{Color.CYAN}{'='*60}{Color.END}")

if __name__ == "__main__":
    import platform
    if platform.system() == "Windows":
        os.system("")
    main()