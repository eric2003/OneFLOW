# test_all_boundaries_fixed.py
import numpy as np
from config import CfdConfig
from mesh import Mesh
from domain import Domain
from solution import Solution
from boundary import BoundaryConditionFactory

print("=== 测试所有边界条件 (修复版) ===")

config = CfdConfig()
mesh = Mesh()
domain = Domain(config, mesh)
solution = Solution(config, domain)

# 测试每种边界条件
boundary_types = ['periodic', 'dirichlet', 'neumann']

for bc_type in boundary_types:
    print(f"\n测试 {bc_type.upper()} 边界:")
    
    # 创建新的配置对象，避免属性污染
    config = CfdConfig()
    config.boundary_type = bc_type
    
    # ✅ 确保设置正确的属性（不是字典键）
    if bc_type == 'dirichlet':
        config.left_boundary_value = 1.0  # 直接设置属性
        config.right_boundary_value = 2.0  # 直接设置属性
    
    mesh = Mesh()
    domain = Domain(config, mesh)
    
    # 创建模拟CFD对象
    class MockCFD:
        def __init__(self, config, domain):
            self.config = config
            self.domain = domain
    
    cfd = MockCFD(config, domain)
    
    # 创建边界条件
    try:
        bc = BoundaryConditionFactory.create(cfd)
        print(f"  ✅ 创建: {type(bc).__name__}")
        
        # 应用边界条件
        u = np.ones(domain.ntcells) * 5.0
        bc.apply(u)
        print(f"  ✅ 应用成功")
        
        # 简单验证
        if bc_type == 'periodic':
            print(f"    左ghost层: {u[:domain.nghosts]}")
            print(f"    右ghost层: {u[-domain.nghosts:]}")
        elif bc_type == 'dirichlet':
            print(f"    左边界值: {u[:domain.nghosts]} (应接近1.0)")
            print(f"    右边界值: {u[-domain.nghosts:]} (应接近2.0)")
        elif bc_type == 'neumann':
            print(f"    左边界梯度: {u[domain.nghosts-1]} = {u[domain.nghosts]}")
            print(f"    右边界梯度: {u[-domain.nghosts]} = {u[-domain.nghosts-1]}")
            
    except Exception as e:
        print(f"  ❌ 失败: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

print("\n=== 测试完成 ===")