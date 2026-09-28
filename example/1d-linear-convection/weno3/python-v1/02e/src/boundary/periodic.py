# src/boundary/periodic.py
from .base import BoundaryCondition
from core.registry import register_component

@register_component('boundary', 'periodic')
class PeriodicBoundary(BoundaryCondition):
    """周期边界条件（1D专用）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        # 左ghost层 = 右物理层
        for ig in range(nghosts):
            u[ist - 1 - ig] = u[ied - 1 - ig]
        
        # 右ghost层 = 左物理层
        for ig in range(nghosts):
            u[ied + ig] = u[ist + ig]