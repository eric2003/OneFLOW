# flux/rusanov.py
"""
Rusanov（Lax-Friedrichs）通量计算器
"""

from cfd_registry import register_component
from .base import InviscidFluxCalculator

@register_component('flux', 'rusanov')
class RusanovFluxCalculator(InviscidFluxCalculator):
    """Rusanov（Lax-Friedrichs）通量"""
    def compute(self, q_face_left, q_face_right, flux):
        for i in range(self.mesh.nnodes):
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            c_L = self.wave_speed
            c_R = self.wave_speed
            F_L = c_L * u_L  # Flux from left state
            F_R = c_R * u_R  # Flux from right state
            Smax = max(abs(c_L), abs(c_R))  # Maximum wave speed
            flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L)