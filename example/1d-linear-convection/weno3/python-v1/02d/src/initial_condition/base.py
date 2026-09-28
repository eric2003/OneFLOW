# src/initial_conditions/base.py
from abc import ABC, abstractmethod
import numpy as np

class InitialCondition(ABC):
    """初始条件抽象基类"""
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
        """内部辅助方法：将值应用到物理区域"""
        domain = solution.domain
        for i in range(domain.ist, domain.ied):
            j = i - domain.ist
            solution.u[i] = values[j]