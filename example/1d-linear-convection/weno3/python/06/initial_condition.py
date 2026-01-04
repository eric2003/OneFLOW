# initial_condition.py

"""
初始条件模块 (已集成注册系统)
使用装饰器 @register_component 自动注册
"""

import numpy as np
from abc import ABC, abstractmethod
from cfd_registry import ComponentRegistry, register_component


# ---------------------- 1. 初始条件抽象基类 ----------------------
class InitialCondition(ABC):
    """初始条件基类"""
    def __init__(self, config):
        self.config = config

    @abstractmethod
    def apply(self, solution):
        """将初始条件应用到 solution 的内部区域"""
        pass

    @abstractmethod
    def evaluate_at(self, x):
        """纯数学函数：给定 x，返回 u(x)，不涉及网格或边界"""
        pass

    def _apply_to_interior(self, solution, values):
        domain = solution.domain
        for i in range(domain.ist, domain.ied):
            j = i - domain.ist
            solution.u[i] = values[j]


# ---------------------- 2. 具体初始条件实现（使用装饰器注册） ----------------------

@register_component('initial_condition', 'step')
class StepFunctionIC(InitialCondition):
    def evaluate_at(self, x):
        u0 = np.ones_like(x)
        mask = (x >= 0.5) & (x <= 1.0)
        u0[mask] = 2.0
        return u0

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)

@register_component('initial_condition', 'sin')
class SineWaveIC(InitialCondition):
    def evaluate_at(self, x):
        L = self.config.get("domain_length", 2.0)
        return np.sin(2 * np.pi * x / L)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)

@register_component('initial_condition', 'gaussian')
class GaussianPulseIC(InitialCondition):
    def evaluate_at(self, x):
        center = self.config.get("pulse_center", 0.5)
        width = self.config.get("pulse_width", 0.1)
        return np.exp(-((x - center) / width) ** 2)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)

class InitialConditionFactory:
    """初始条件工厂"""
    
    @classmethod
    def create(cls, ic_type: str, config) -> 'InitialCondition':
        """
        创建初始条件实例
        
        Args:
            ic_type: 初始条件类型（如'step', 'sin'等）
            config: 配置对象
        
        Returns:
            初始条件实例
        """
        from factories.base_factory import BaseFactory
        return BaseFactory.create_component('initial_condition', ic_type, config)