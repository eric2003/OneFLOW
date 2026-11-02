# initial_condition.py
import numpy as np
from abc import ABC, abstractmethod

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


# ---------------------- 2. 具体初始条件实现 ----------------------
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


class SineWaveIC(InitialCondition):
    def evaluate_at(self, x):
        L = self.config.get("domain_length", 2.0)
        return np.sin(2 * np.pi * x / L)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)


class GaussianPulseIC(InitialCondition):
    def evaluate_at(self, x):
        center = self.config.get("pulse_center", 0.5)
        width = self.config.get("pulse_width", 0.1)
        return np.exp(-((x - center) / width) ** 2)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)


# ---------------------- 3. 初始条件工厂 ----------------------
class InitialConditionFactory:
    _registry = {
        'step': StepFunctionIC,
        'sin': SineWaveIC,
        'gaussian': GaussianPulseIC,
    }

    @classmethod
    def create(cls, ic_type, config):
        if ic_type not in cls._registry:
            raise ValueError(f"未知的初始条件类型: {ic_type}（支持: {list(cls._registry.keys())}）")
        return cls._registry[ic_type](config)

    @classmethod
    def register(cls, name, ic_class):
        cls._registry[name] = ic_class