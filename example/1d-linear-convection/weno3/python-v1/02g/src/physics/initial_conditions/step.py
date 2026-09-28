# src/initial_conditions/step.py
import numpy as np
from .base import InitialCondition
from core.registry import register_component

@register_component('initial_condition', 'step')
class StepFunctionIC(InitialCondition):
    def evaluate_at(self, x):
        u0 = np.ones_like(x)
        mask = (x >= 0.5) & (x <= 1.0)
        u0[mask] = 2.0
        return u0

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)