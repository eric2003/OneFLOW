# src/physics/problems/linear_advection.py
"""
线性对流问题
"""

from .base import Problem
from physics.equations.linear_advection import LinearAdvectionEquation
from physics.systems.linear_advection_system import LinearAdvectionSystem
from physics.initial_conditions.factory import InitialConditionFactory
from core.registry import register_component

@register_component('problem', 'linear_advection')
class LinearAdvectionProblem(Problem):
    """
    线性对流问题：u_t + c u_x = 0
    使用周期边界 + 任意初始条件
    """

    def __init__(self, config):
        super().__init__(config)
        
    def initial_condition(self):
        """返回初始条件实例"""
        if self._initial_condition is None:
            from physics.initial_conditions.factory import InitialConditionFactory
            self._initial_condition = InitialConditionFactory.create(self.config)
        return self._initial_condition

    def physical_system(self):
        """返回物理系统实例"""
        if self._physical_system is None:
            # 创建方程
            equation = LinearAdvectionEquation(wave_speed=self.config.wave_speed)
            # 创建系统
            self._physical_system = LinearAdvectionSystem(equation)
        return self._physical_system