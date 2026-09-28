# src/initial_conditions/gaussian.py
import numpy as np
from .base import InitialCondition
from core.registry import register_component

@register_component('initial_condition', 'gaussian')
class GaussianPulseIC(InitialCondition):
    def evaluate_at(self, x):
        center = getattr(self.config, "pulse_center", 0.5)
        width = getattr(self.config, "pulse_width", 0.1)
        return np.exp(-((x - center) / width) ** 2)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)