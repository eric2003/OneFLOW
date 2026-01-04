# python/gen_ic_test_data.py
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from config import CfdConfig
from mesh import Mesh
from domain import Domain
from solution import Solution

# 固定 mesh
mesh = Mesh()
config = CfdConfig()

# 测试三种 IC
for ic_type in ["step", "sin", "gaussian"]:
    config.ic_type = ic_type
    domain = Domain(config, mesh)
    sol = Solution(config, domain)
    
    u_full = sol.u.copy()  # 包含 ghost
    u_interior = sol.u[domain.ist:domain.ied].copy()
    
    np.save(f"u_{ic_type}_full_py.npy", u_full)
    np.save(f"u_{ic_type}_interior_py.npy", u_interior)

print("✅ 初始条件测试数据已生成")