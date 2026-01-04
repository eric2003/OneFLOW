# solver.py

from mesh import Mesh
from domain import Domain
from solution import Solution
from config import CfdConfig
from result import ResultAssembler
from problem import Problem
from initial_condition import InitialConditionFactory
from boundary.factory import BoundaryConditionFactory

from numerics.time_integration.factory import TimeIntegratorFactory
from numerics.residual import ResidualCalculator

# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, problem: Problem, mesh):
        self.problem = problem
        self.config = problem.config
        self.domain = Domain(problem.config, mesh)
        self.solution = Solution(problem.config, self.domain)
        
        self.initial_condition = InitialConditionFactory.create(self.config)
        self.boundary_condition = BoundaryConditionFactory.create(self)
        self.residual_calculator = ResidualCalculator(self)
        self.integrator = TimeIntegratorFactory.create(self)
        
    # ==================== 公共接口（供 TimeIntegrator 调用） ====================
    def compute_residual(self):
        """计算物理残差（封装重建→通量→散度）"""
        self.residual_calculator.compute()

    def apply_boundary(self):
        """应用边界条件"""
        self.boundary_condition.apply(self.solution.u)

    # ==================== 初始化 ====================
    def initialize(self):
        """初始化全场：先 IC，再 BC，最后同步 old field"""
        self.initial_condition.apply(self.solution)
        # 应用边界条件到初始场
        self.apply_boundary()
        # 同步 old field
        self.solution.update_old_field()        
        

    def run(self):
        self.initialize()
        
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
        
        self.result = ResultAssembler.assemble(self)
        return self.result["numerical"]
        