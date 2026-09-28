# reconstructor/weno3.py
import numpy as np
from .base import Reconstructor

class Weno3Reconstructor(Reconstructor):
    def reconstruct(self, q, cfd):
        domain = cfd.domain
        solution = cfd.solution
        self._reconstruct_left_interfaces(domain, q, solution.q_face_left)
        self._reconstruct_right_interfaces(domain, q, solution.q_face_right)
        
    def _reconstruct_left_interfaces(self, domain, u, qL):
        """在每个 i+1/2 界面，计算左单元贡献的 qL (即 u_{i+1/2}^-)"""
        for i in range(domain.ist - 1, domain.ied):
            j = i - (domain.ist - 1)
            v1, v2, v3 = u[i-1], u[i], u[i+1]
            qL[j] = self._reconstruct_from_right_biased_stencil(v3, v2, v1)
            
    def _reconstruct_right_interfaces(self, domain, u, qR):
        """在每个 i+1/2 界面，计算右单元贡献的 qR (即 u_{i+1/2}^+)"""
        for i in range(domain.ist, domain.ied + 1):
            j = i - domain.ist
            v1, v2, v3 = u[i-1], u[i], u[i+1]
            qR[j] = self._reconstruct_from_left_biased_stencil(v3, v2, v1)

    def _reconstruct_from_left_biased_stencil(self, v1, v2, v3):
        """使用左偏 stencil (v1,v2,v3) 重建界面值（对应 u_{i+1/2}^+）"""
        eps = 1e-6
        beta0 = (v3 - v2)**2      # smoothness indicator for stencil [v2, v3]
        beta1 = (v2 - v1)**2      # smoothness indicator for stencil [v1, v2]
        
        d0, d1 = 2/3, 1/3         # optimal linear weights (for right value)
        alpha0 = d0 / (eps + beta0)**2
        alpha1 = d1 / (eps + beta1)**2
        w0 = alpha0 / (alpha0 + alpha1)
        w1 = alpha1 / (alpha0 + alpha1)
        
        q0 = 1.5 * v2 - 0.5 * v3  # reconstruction from [v2, v3]
        q1 = 0.5 * v1 + 0.5 * v2  # reconstruction from [v1, v2]
        return w0 * q0 + w1 * q1

    def _reconstruct_from_right_biased_stencil(self, v1, v2, v3):
        """使用右偏 stencil (v1,v2,v3) 重建界面值（对应 u_{i+1/2}^-）"""
        eps = 1e-6
        beta0 = (v3 - v2)**2
        beta1 = (v2 - v1)**2
        
        d0, d1 = 1/3, 2/3         # optimal linear weights (for left value)
        alpha0 = d0 / (eps + beta0)**2
        alpha1 = d1 / (eps + beta1)**2
        w0 = alpha0 / (alpha0 + alpha1)
        w1 = alpha1 / (alpha0 + alpha1)
        
        q0 = 1.5 * v2 - 0.5 * v3  # from [v2, v3]
        q1 = 0.5 * v1 + 0.5 * v2 # from [v1, v2]
        return w0 * q0 + w1 * q1
