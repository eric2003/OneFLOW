# src/physics/systems/base.py
"""
物理系统基类
将方程与解析解等逻辑组合
"""

from abc import ABC, abstractmethod

class PhysicalSystem(ABC):
    """物理系统抽象基类"""
    
    def __init__(self, equation):
        self.equation = equation
    
    @abstractmethod
    def exact_solution(self, x, t, initial_condition):
        """
        计算给定初始条件的解析解
        :param x: 空间坐标数组
        :param t: 时间
        :param initial_condition: 初始条件对象
        :return: 解析解数组
        """
        pass
    
    @property
    def wave_speed(self):
        """获取波速（委托给方程）"""
        return self.equation.wave_speed()
    
    def flux(self, u):
        """计算通量（委托给方程）"""
        return self.equation.flux(u)