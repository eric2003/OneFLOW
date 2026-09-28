# src/boundary/dirichlet.py
from .base import BoundaryCondition
from core.registry import register_component

@register_component('boundary', 'dirichlet')
class DirichletBoundary(BoundaryCondition):
    """Dirichlet（固定值）边界条件（如进口固定速度、固壁零速度）"""
    def apply(self, u):
        nghosts = self.domain.nghosts
        ist = self.domain.ist
        ied = self.domain.ied
        
        left_value = getattr(self.config, "left_boundary_value", 1.0)
        
        # 左边界（进口）固定值
        for ig in range(nghosts):
            u[ist - 1 - ig] = left_value
        
        right_value = getattr(self.config, "right_boundary_value", 2.0)
        
        # 右边界（出口）固定值
        for ig in range(nghosts):
            u[ied + ig] = right_value
        
        # 调试信息
        if hasattr(self.config, 'debug') and self.config.debug:
            print(f"  应用Dirichlet边界: 左值={left_value}, 右值={right_value}")