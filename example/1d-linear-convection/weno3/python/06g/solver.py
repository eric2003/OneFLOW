from boundary import BoundaryConditionFactory
from initial_condition import InitialConditionFactory
from time_integration import TimeIntegrator,TimeIntegratorFactory
from mesh import Mesh
from domain import Domain
from solution import Solution
from config import CfdConfig
from result import ResultAssembler

# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, config, mesh):
        self.config = config
        self.domain = Domain(config, mesh)
        self.solution = Solution(config, self.domain)
        self.integrator = TimeIntegratorFactory.create(self)
        self.boundary_condition = BoundaryConditionFactory.create(self)

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
        
        self.result = ResultAssembler.assemble(self)
        return self.result["numerical"]        
        