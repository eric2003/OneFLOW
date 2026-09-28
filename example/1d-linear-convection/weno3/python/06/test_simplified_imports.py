# test_simplified_imports.py
"""
测试简化后的导入系统
"""

print("=== 测试简化后的导入系统 ===")

# 测试导入所有模块
modules = [
    'boundary',
    'flux', 
    'initial_condition',
    'time_integration',
]

print("\n1. 导入所有模块:")
for module in modules:
    try:
        __import__(module)
        print(f"  ✅ {module}")
    except ImportError as e:
        print(f"  ❌ {module}: {e}")

# 测试注册系统
print("\n2. 测试注册系统:")
from cfd_registry import ComponentRegistry

print("所有注册的组件:")
all_components = ComponentRegistry.list_all()
for category, names in all_components.items():
    print(f"  {category}: {names}")

print(f"\n总计: {ComponentRegistry.get_count()} 个组件")

# 测试关闭verbose模式
print("\n3. 测试关闭verbose模式:")
ComponentRegistry.set_verbose(False)

# 尝试重新导入一个模块（应该不会打印注册信息）
print("重新导入boundary模块（应该静默）:")
import importlib
import boundary
importlib.reload(boundary)

print("\n✅ 简化导入系统测试完成!")