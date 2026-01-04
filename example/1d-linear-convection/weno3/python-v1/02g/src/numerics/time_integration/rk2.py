# time_integration/rk2.py

from .base import TimeIntegrator
from core.registry import register_component

@register_component('integrator', 'rk2')
class RK2Integrator(TimeIntegrator):
    """2阶Heun方法（RK2）"""
    def step(self, dt):
        # 预测步
        self.compute_residual()
        u_pred = self.solution.u.copy()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u_pred[i] += dt * self.solution.res[j]
        self.solution.u[:] = u_pred
        self.apply_boundary()

        # 校正步
        self.compute_residual()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] = (
                0.5 * self.solution.un[i] +
                0.5 * self.solution.u[i] +
                0.5 * dt * self.solution.res[j]
            )
        self.apply_boundary()
        self.solution.update_old_field()