# flux/rusanov.py
"""
Rusanov（Lax-Friedrichs）通量计算器
"""

from core.registry import register_component
from .base import InviscidFluxCalculator

@register_component('flux', 'rusanov')
class RusanovFluxCalculator(InviscidFluxCalculator):
    """Rusanov（Lax-Friedrichs）通量，使用 Equation 解耦物理参数"""

    def compute(self, q_face_left, q_face_right, flux):
        # 从 cfd 获取 equation（与 Julia 对齐）
        eq = self.cfd.problem.physical_system().equation
        for i in range(self.mesh.nnodes):
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            c_L = eq.wave_speed()
            c_R = eq.wave_speed()
            # 通过 equation 计算通量和波速
            F_L = eq.flux(u_L)
            F_R = eq.flux(u_R)
            Smax = max(abs(c_L), abs(c_R))  # Maximum wave speed
            flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L)            