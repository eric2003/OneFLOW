#!/usr/bin/env python3
# scripts/test_integrated.py
"""
测试集成求解器
"""
import subprocess
import sys
import os
from pathlib import Path

def run_test():
    """运行集成测试"""
    # 构建项目
    print("构建项目...")
    result = subprocess.run(
        ["python", "scripts/build.py", "--no-tests", "--clean"],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print("构建失败:")
        print(result.stderr)
        return False
    
    # 运行集成示例
    print("\n运行集成示例...")
    exe_path = Path("build/bin/Debug/run_eno_weno_integrated.exe")
    
    if not exe_path.exists():
        print(f"可执行文件不存在: {exe_path}")
        return False
    
    result = subprocess.run(
        [str(exe_path)],
        capture_output=True,
        text=True,
        encoding='utf-8'
    )
    
    print(result.stdout)
    if result.stderr:
        print("错误输出:")
        print(result.stderr)
    
    return result.returncode == 0

if __name__ == "__main__":
    success = run_test()
    sys.exit(0 if success else 1)