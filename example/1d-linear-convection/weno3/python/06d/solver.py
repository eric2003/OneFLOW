import numpy as np
import matplotlib.pyplot as plt
import inflect
from abc import ABC, abstractmethod

from boundary import BoundaryConditionFactory
from initial_condition import InitialConditionFactory
from time_integration import TimeIntegrator,TimeIntegratorFactory
from flux import InviscidFluxCalculator  # 仅用于类型提示（可选）
from mesh import Mesh
from domain import Domain
from solution import Solution
from config import CfdConfig
from residual import ResidualCalculator
from reconstructor import ReconstructorFactory
from factories.base_factory import BaseFactory

# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, config, mesh):
        self.config = config
        self.domain = Domain(config, mesh)
        self.solution = Solution(config, self.domain)
        self.reconstructor = ReconstructorFactory.create(config, self.domain)
        self.residual_calculator = ResidualCalculator(self)
        self.integrator = TimeIntegratorFactory.create(self)
        self.boundary_condition = BoundaryConditionFactory.create(self)
   
    
    def exact_solution(self):
        """通用对流问题的解析解：u(x, T) = u0(x - c*T)，周期边界"""
        x = self.domain.mesh.xcc
        T = self.config.final_time
        c = self.config.wave_speed
        L = self.domain.mesh.L
 
        # 周期平移：确保在 [x0, x0 + L) 内
        x_shifted = (x - c * T + L) % L

        # 获取 IC 实例并评估
        ic = InitialConditionFactory.create(self.config)
        return ic.evaluate_at(x_shifted)
    
    def run(self):
        # 应用初始边界条件并同步 old field
        self.boundary_condition.apply(self.solution.u)
        self.solution.update_old_field()       
        
        t = 0.0
        dt_old = self.config.dt
        dt = dt_old
        
        domain = self.domain
        solution = self.solution
        config = self.config
        
        while t < config.final_time:
            if t + dt > config.final_time:
                dt = config.final_time - t
            config.dt = dt  # temporary adjustment for last step
            self.integrator.step(dt)
            t += dt
        config.dt = dt_old
        
        # 整理标准化结果
        u_numerical = self.solution.u[self.domain.ist:self.domain.ied].copy()
        self.result = {
            "x": domain.mesh.xcc,
            "numerical": u_numerical,
            "analytical": self.exact_solution(),
            "config": {
                "scheme": self.config.recon_scheme,
                "order": self.config.spatial_order,
                "rk_order": self.config.rk_order,
                "final_time": self.config.final_time
            }
        }

        return u_numerical
