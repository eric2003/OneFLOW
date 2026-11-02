# src/physics/systems/linear_advection_system.py
"""
线性对流物理系统
"""

from .base import PhysicalSystem
from core.registry import register_component
import numpy as np

@register_component('system', 'linear_advection')
class LinearAdvectionSystem(PhysicalSystem):
    """线性对流系统"""
    
    def exact_solution(self, x, t, initial_condition):
        """
        计算线性对流的解析解
        :param x: 空间坐标数组
        :param t: 时间
        :param initial_condition: 初始条件对象
        :return: 解析解数组
        """
        if len(x) == 0:
            return np.zeros((1, 0))
        
        L = x[-1] - x[0]  # 假设均匀周期网格
        c = self.equation.wave_speed()
        x_shifted = (x - c * t + L) % L
        
        # 调用初始条件的 evaluate_at 方法
        u0 = initial_condition.evaluate_at(x_shifted)
        
        # 确保返回形状正确
        if u0.ndim == 1:
            u0 = u0[np.newaxis, :]  # (1, ncells)
        
        return u0