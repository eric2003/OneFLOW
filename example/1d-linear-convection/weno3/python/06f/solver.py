from boundary import BoundaryConditionFactory
from initial_condition import InitialConditionFactory
from time_integration import TimeIntegrator,TimeIntegratorFactory
from mesh import Mesh
from domain import Domain
from solution import Solution
from config import CfdConfig

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
        
        self.result = self._assemble_result()
        return self.result["numerical"]
        
    def _assemble_result(self):
        """组装标准化求解结果字典"""
        u_numerical = self.solution.u[self.domain.ist:self.domain.ied].copy()
        ic = InitialConditionFactory.create(self.config)
        analytical = ic.exact_solution(self)
        return {
            "x": self.domain.mesh.xcc,
            "numerical": u_numerical,
            "analytical": analytical,
            "config": {
                "scheme": self.config.recon_scheme,
                "order": self.config.spatial_order,
                "rk_order": self.config.rk_order,
                "final_time": self.config.final_time
            }
        }