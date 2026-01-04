# src/physics/equations/euler.py
"""
欧拉方程（占位文件，未来实现）
"""
from .base import Equation

class EulerEquation(Equation):
    """欧拉方程：质量、动量、能量守恒"""
    def __init__(self, gamma=1.4):
        self.gamma = gamma
    
    @property
    def num_equations(self) -> int:
        return 3
    
    def flux(self, u):
        raise NotImplementedError("欧拉方程待实现")
    
    def max_wave_speed(self, u):
        raise NotImplementedError("欧拉方程待实现")
        
    def wave_speed(self):
        raise NotImplementedError("欧拉方程待实现")