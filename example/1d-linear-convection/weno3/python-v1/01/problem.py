# problem.py
"""
问题定义模块：每个 Problem 子类代表一个完整测试用例
包含初始条件、解析解（可选）、物理方程等
"""

from abc import ABC, abstractmethod
from initial_condition import InitialConditionFactory


class Problem(ABC):
    """
    抽象问题基类
    每个具体问题（如线性对流、Sod 激波管）应继承此类
    """
    def __init__(self, config):
        self.config = config

    @abstractmethod
    def initial_condition(self):
        """
        返回 InitialCondition 实例
        """
        pass

    def exact_solution(self, cfd):
        """
        可选：返回解析解（数值数组）
        若无解析解，可抛出 NotImplementedError 或返回 None
        """
        x = cfd.domain.mesh.xcc
        raise NotImplementedError(
            f"Problem '{self.__class__.__name__}' does not provide an exact solution."
        )

    @property
    def name(self):
        return self.__class__.__name__