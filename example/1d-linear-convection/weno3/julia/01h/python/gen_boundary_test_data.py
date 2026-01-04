# python/gen_boundary_test_data.py
import numpy as np
from config import CfdConfig
from mesh import Mesh
from domain import Domain

# 固定测试配置
config = CfdConfig()
config.with_boundary("dirichlet", left_value=0.5, right_value=1.5)
config.debug = True

mesh = Mesh()
domain = Domain(config, mesh)

# 构造 mock CFD 对象（仅含 config + domain）
class MockCfd:
    def __init__(self, config, domain):
        self.config = config
        self.domain = domain

# 测试用 u：0,1,2,...,N-1
u_input = np.arange(domain.ntcells, dtype=np.float64)
np.save("u_input.npy", u_input)

# 测试每种边界
from boundary import PeriodicBoundary, DirichletBoundary, NeumannBoundary

for bc_name, bc_class in [("periodic", PeriodicBoundary), ("dirichlet", DirichletBoundary), ("neumann", NeumannBoundary)]:
    u = u_input.copy()
    cfd_mock = MockCfd(config, domain)
    bc = bc_class(cfd_mock)
    bc.apply(u)
    np.save(f"u_{bc_name}_py.npy", u)

print("✅ 测试数据已生成：u_*.npy")