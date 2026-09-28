# src/initial_conditions/sine.py
import numpy as np
from .base import InitialCondition
from core.registry import register_component

@register_component('initial_condition', 'sin')
class SineWaveIC(InitialCondition):
    def evaluate_at(self, x):
        L = getattr(self.config, "domain_length", 2.0)
        return np.sin(2 * np.pi * x / L)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)