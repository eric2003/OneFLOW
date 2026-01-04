# flux/engquist_osher.py
"""
Engquist-Osher 通量计算器（线性对流专用）
"""

from core.registry import register_component
from .base import InviscidFluxCalculator

@register_component('flux', 'engquist-osher')
class EngquistOsherFluxCalculator(InviscidFluxCalculator):
    """Engquist-Osher通量（线性对流专用）"""
    def compute(self, q_face_left, q_face_right, flux):
        for i in range(self.mesh.nnodes):
            c = self.wave_speed
            cp = 0.5 * (c + abs(c))
            cm = 0.5 * (c - abs(c))
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            flux[i] = cp * u_L + cm * u_R