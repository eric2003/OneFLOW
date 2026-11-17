# time_integration/base.py

from abc import ABC, abstractmethod

class TimeIntegrator(ABC):
    """时间推进器抽象基类：定义一维CFD时间推进的核心接口"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.config = cfd.config
        self.domain = cfd.domain
        self.solution = cfd.solution

    @abstractmethod
    def step(self, dt):
        pass

    def compute_residual(self):
        """计算残差（委托给 Cfd）"""
        self.cfd.compute_residual()

    def apply_boundary(self):
        """应用边界条件（委托给 Cfd）"""
        self.cfd.apply_boundary()

    def map_idx(self, i):
        """物理网格索引 → 残差数组索引"""
        return i - self.domain.ist