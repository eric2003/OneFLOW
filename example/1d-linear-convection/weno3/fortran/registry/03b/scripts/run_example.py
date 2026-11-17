#!/usr/bin/env python3
# scripts/run_example.py
"""
运行ENO/WENO示例程序
"""

import os
import sys
import subprocess
from pathlib import Path

# 添加当前目录到路径
sys.path.insert(0, str(Path(__file__).parent))

try:
    from build import BuildSystem
except ImportError:
    print("错误: 找不到build.py，请确保在scripts目录中运行")
    sys.exit(1)

def run_example():
    """运行示例程序"""
    builder = BuildSystem()
    
    print("\n" + "="*70)
    print("  运行ENO/WENO对比示例程序")
    print("="*70 + "\n")
    
    # 检查是否已构建
    exe_path = builder.build_dir / "bin" / "Debug" / "example_eno_weno_comparison.exe"
    
    if not exe_path.exists():
        print("示例程序未构建，先构建项目...")
        print("-"*50)
        
        # 使用简化的构建
        result = subprocess.run(
            ["python", "build.py", "--no-tests", "--clean"],
            cwd=builder.project_root / "scripts",
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print("构建失败:")
            print(result.stderr)
            return False
    
    # 运行示例程序
    print("运行示例程序...")
    print("-"*50)
    
    try:
        result = subprocess.run(
            [str(exe_path)],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace'
        )
        
        print(result.stdout)
        
        if result.stderr:
            print("标准错误输出:")
            print(result.stderr)
        
        return result.returncode == 0
        
    except Exception as e:
        print(f"运行示例程序失败: {e}")
        return False

def main():
    """主函数"""
    success = run_example()
    
    if success:
        print("\n" + "="*70)
        print("  ✓ 示例程序运行成功")
        print("="*70)
        return 0
    else:
        print("\n" + "="*70)
        print("  ✗ 示例程序运行失败")
        print("="*70)
        return 1

if __name__ == "__main__":
    sys.exit(main())