# src/physics/problems/base.py
"""
问题定义抽象基类
每个具体问题（如线性对流、Sod 激波管）应继承此类
"""

from abc import ABC, abstractmethod

class Problem(ABC):
    """
    抽象问题基类
    """
    def __init__(self, config):
        self.config = config
        self._initial_condition = None
        self._physical_system = None

    @abstractmethod
    def initial_condition(self):
        """返回 InitialCondition 实例"""
        pass
    
    @abstractmethod
    def physical_system(self):
        """返回 PhysicalSystem 实例"""
        pass

    def exact_solution(self, cfd):
        """
        计算解析解（委托给物理系统）
        """
        x = cfd.domain.mesh.xcc
        t = cfd.config.final_time
        ic = self.initial_condition()
        system = self.physical_system()
        
        u_exact = system.exact_solution(x, t, ic)
        return u_exact[0]  # 返回标量数组 (ncells,) 以兼容现有绘图逻辑

    @property
    def name(self):
        return self.__class__.__name__