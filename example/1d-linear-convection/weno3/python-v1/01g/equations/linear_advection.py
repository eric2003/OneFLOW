# equations/linear_advection.py

import numpy as np
from .base import Equation

class LinearAdvectionEquation(Equation):
    """
    线性对流方程: u_t + c * u_x = 0
    支持标量（num_equations=1）
    """

    def __init__(self, wave_speed: float):
        self.c = wave_speed

    @property
    def num_equations(self) -> int:
        return 1

    def flux(self, u):
        """f(u) = c * u"""
        return self.c * u

    def max_wave_speed(self, u):
        """最大波速 = |c|"""
        return abs(self.c)
        
    def wave_speed(self):
        return self.c
        

    def exact_solution(self, x, t, initial_condition_func):
        """
        解析解: u(x, t) = u0(x - c * t)
        支持周期边界
        """
        if len(x) == 0:
            return np.zeros((1, 0))
        
        L = x[-1] - x[0]  # 假设均匀周期网格
        x_shifted = (x - self.c * t + L) % L
        u0 = initial_condition_func(x_shifted)  # (1, ncells)
        return u0