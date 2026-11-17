# time_integration/rk1.py

from .base import TimeIntegrator
from core.registry import register_component

@register_component('integrator', 'rk1')
class RK1Integrator(TimeIntegrator):
    """1阶显式欧拉（RK1）"""
    def step(self, dt):
        self.compute_residual()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] += dt * self.solution.res[j]
        self.apply_boundary()
        self.solution.update_old_field()