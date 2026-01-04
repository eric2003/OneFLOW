# time_integration/rk3.py

from .base import TimeIntegrator
from cfd_registry import register_component

@register_component('integrator', 'rk3')
class RK3Integrator(TimeIntegrator):
    """3阶SSPRK3（强稳定保号RK3）"""
    def step(self, dt):
        # Stage 1
        self.compute_residual()
        u1 = self.solution.u.copy()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u1[i] += dt * self.solution.res[j]
        self.solution.u[:] = u1
        self.apply_boundary()

        # Stage 2
        self.compute_residual()
        u2 = self.solution.u.copy()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u2[i] = (
                0.75 * self.solution.un[i] +
                0.25 * self.solution.u[i] +
                0.25 * dt * self.solution.res[j]
            )
        self.solution.u[:] = u2
        self.apply_boundary()

        # Stage 3
        self.compute_residual()
        c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] = (
                c1 * self.solution.un[i] +
                c2 * self.solution.u[i] +
                c3 * dt * self.solution.res[j]
            )
        self.apply_boundary()
        self.solution.update_old_field()