# test_registry_complete.py
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("=== 完整的注册系统测试 ===")

# 导入所有模块以触发注册
print("\n1. 导入所有CFD模块（触发注册）:")
modules = [
    'boundary',
    'flux', 
    'initial_condition',
    'time_integration',
    # reconstructor 已经在导入时注册
]

for module in modules:
    try:
        __import__(module)
        print(f"   ✅ {module}")
    except ImportError as e:
        print(f"   ❌ {module}: {e}")

# 检查注册表
print("\n2. 检查注册表内容:")

from registry import ComponentRegistry

all_components = ComponentRegistry.list_all()
total_components = 0

for category, components in all_components.items():
    print(f"\n  {category.upper()} ({len(components)}种):")
    for name in sorted(components):
        try:
            cls = ComponentRegistry.get(category, name)
            print(f"    - {name}: {cls.__name__}")
            total_components += 1
        except:
            print(f"    - {name}: <无法获取类>")

print(f"\n总计: {total_components} 个组件已注册")
print("\n3. 测试创建每个类别的组件:")

# 测试数据
class MockConfig:
    boundary_type = 'periodic'
    flux_type = 'rusanov'
    recon_scheme = 'eno'
    rk_order = 2
    ic_type = 'step'
    spatial_order = 3

class MockCFD:
    def __init__(self):
        self.config = MockConfig()
        self.domain = type('Domain', (), {'ntcells': 50, 'mesh': type('Mesh', (), {'nnodes': 10})()})()

# 测试每个类别
test_cases = [
    ('boundary', 'periodic', MockCFD()),
    ('flux', 'rusanov', MockCFD()),
    ('initial_condition', 'step', MockConfig()),
    ('integrator', 'rk2', MockCFD()),
]

for category, name, arg in test_cases:
    try:
        instance = ComponentRegistry.create(category, name, arg)
        print(f"  ✅ {category}.{name}: 创建成功 -> {type(instance).__name__}")
    except Exception as e:
        print(f"  ❌ {category}.{name}: 创建失败 -> {type(e).__name__}: {e}")

print("\n✅ 注册系统测试完成!")