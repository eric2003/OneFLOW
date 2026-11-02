# src/boundary/neumann.py
from .base import BoundaryCondition
from core.registry import register_component

@register_component('boundary', 'neumann')
class NeumannBoundary(BoundaryCondition):
    """Neumann（零梯度）边界条件（如出口无梯度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左边界零梯度
        for ig in range(nghosts):
            u[ist - 1 - ig] = u[ist + ig]
        
        # 右边界零梯度
        for ig in range(nghosts):
            u[ied + ig] = u[ied - 1 - ig]
        
        # 调试信息
        if hasattr(self.config, 'debug') and self.config.debug:
            print(f"  应用Neumann边界: 零梯度")