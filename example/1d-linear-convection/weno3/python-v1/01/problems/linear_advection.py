# problems/linear_advection.py
"""
线性对流问题：支持任意初始条件 + 周期平移解析解
"""

import numpy as np
from problem import Problem
from initial_condition import InitialConditionFactory


class LinearAdvectionProblem(Problem):
    """
    线性对流问题 u_t + c u_x = 0
    解析解：u(x, t) = u0(x - c * t)
    """

    def initial_condition(self):
        return InitialConditionFactory.create(self.config)

    def exact_solution(self, cfd):
        ic = self.initial_condition()
        x = cfd.domain.mesh.xcc
        t = cfd.config.final_time
        c = cfd.config.wave_speed
        L = cfd.domain.mesh.L
        x_shifted = (x - c * t + L) % L
        return ic.evaluate_at(x_shifted)