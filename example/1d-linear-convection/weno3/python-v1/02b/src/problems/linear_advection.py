# problems/linear_advection.py

import numpy as np
from problem import Problem
from equations.linear_advection import LinearAdvectionEquation
from initial_condition import InitialConditionFactory

class LinearAdvectionProblem(Problem):
    """
    线性对流问题：u_t + c u_x = 0
    使用周期边界 + 任意初始条件
    """

    def __init__(self, config):
        super().__init__(config)
        # 持有物理方程实例
        self.equation = LinearAdvectionEquation(wave_speed=config.wave_speed)
        self._ic_cache = None  # 缓存 IC 实例

    def initial_condition(self):
        """返回初始条件实例"""
        if self._ic_cache is None:
            self._ic_cache = InitialConditionFactory.create(self.config)
        return self._ic_cache

    def exact_solution(self, cfd):
        """
        委托给 Equation 计算解析解
        """
        ic = self.initial_condition()
        # 包装 evaluate_at 以返回 (1, ncells)
        def vectorized_ic(x):
            u0_scalar = ic.evaluate_at(x)  # (ncells,)
            return u0_scalar[np.newaxis, :]  # (1, ncells)
        
        x = cfd.domain.mesh.xcc
        t = cfd.config.final_time
        u_exact = self.equation.exact_solution(x, t, vectorized_ic)
        return u_exact[0]  # 返回标量数组 (ncells,) 以兼容现有绘图逻辑