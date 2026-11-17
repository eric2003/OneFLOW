# reconstructor/weno5.py
import numpy as np
from .base import Reconstructor
from core.registry import ComponentRegistry, register_component

@register_component('reconstructor', 'weno5')
class Weno5Reconstructor(Reconstructor):
    def compute_face_values(self, q, cfd):
        domain = cfd.domain
        solution = cfd.solution
        self._reconstruct_left_interfaces(domain, q, solution.q_face_left)
        self._reconstruct_right_interfaces(domain, q, solution.q_face_right)
        
    def _reconstruct_left_interfaces(self, domain, u, qL):
        """在每个 i+1/2 界面，计算左单元贡献的 qL (即 u_{i+1/2}^-)"""
        for i in range(domain.ist - 1, domain.ied):
            j = i - (domain.ist - 1)
            v1, v2, v3 = u[i-2], u[i-1], u[i], u[i+1], u[i+2]
            #qL[j] = self._reconstruct_from_right_biased_stencil(v3, v2, v1)
            qL[j] = self._reconstruct_from_left_biased_stencil(v1, v2, v3, v4, v5)
            
    def _reconstruct_right_interfaces(self, domain, u, qR):
        """在每个 i+1/2 界面，计算右单元贡献的 qR (即 u_{i+1/2}^+)"""
        for i in range(domain.ist, domain.ied + 1):
            j = i - domain.ist
            v1, v2, v3 = u[i-2], u[i-1], u[i], u[i+1], u[i+2]
            #qR[j] = self._reconstruct_from_left_biased_stencil(v3, v2, v1)
            qR[j] = self._reconstruct_from_right_biased_stencil(v1, v2, v3, v4, v5)
            
    def _reconstruct_from_left_biased_stencil(self, v1, v2, v3, v4, v5):
        eps = 1e-6
        beta0 = (13.0/12.0)*(v1 - 2*v2 + v3)**2 + (1.0/4.0)*(v1 - 4*v2 + 3*v3)**2
        beta1 = (13.0/12.0)*(v2 - 2*v3 + v4)**2 + (1.0/4.0)*(v2 - v4)**2
        beta2 = (13.0/12.0)*(v3 - 2*v4 + v5)**2 + (1.0/4.0)*(3*v3 - 4*v4 + v5)**2
        d0 = 1.0/10.0
        d1 = 3.0/5.0
        d2 = 3.0/10.0
        alpha0 = d0 / (eps + beta0)**2
        alpha1 = d1 / (eps + beta1)**2
        alpha2 = d2 / (eps + beta2)**2
        alpha = alpha0 + alpha1 + alpha2
        w0 = alpha0 / alpha
        w1 = alpha1 / alpha
        w2 = alpha2 / alpha
        q0 = 1.0/3.0*v1-7.0/6.0*v2+11.0/6.0*v3  # r=2
        q1 = -1.0/6.0*v2+5.0/6.0*v3+1.0/3.0*v4  # r=1
        q2 = 1.0/3.0*v3+5.0/6.0*v4-1.0/6.0*v5  # r=0
        return w0 * q0 + w1 * q1 + w2 * q2

    def _reconstruct_from_right_biased_stencil(self, v1, v2, v3, v4, v5):
        eps = 1e-6
        beta0 = (13.0/12.0)*(v1 - 2*v2 + v3)**2 + (1.0/4.0)*(v1 - 4*v2 + 3*v3)**2
        beta1 = (13.0/12.0)*(v2 - 2*v3 + v4)**2 + (1.0/4.0)*(v2 - v4)**2
        beta2 = (13.0/12.0)*(v3 - 2*v4 + v5)**2 + (1.0/4.0)*(3*v3 - 4*v4 + v5)**2
        d0 = 3.0/10.0
        d1 = 3.0/5.0
        d2 = 1.0/10.0
        alpha0 = d0 / (eps + beta0)**2
        alpha1 = d1 / (eps + beta1)**2
        alpha2 = d2 / (eps + beta2)**2
        alpha = alpha0 + alpha1 + alpha2
        w0 = alpha0 / alpha
        w1 = alpha1 / alpha
        w2 = alpha2 / alpha
        q0 = -1.0/6.0*v1+5.0/6.0*v2+1.0/3.0*v3  # r=2
        q1 = 1.0/3.0*v2+5.0/6.0*v3-1.0/6.0*v4  # r=1
        q2 = 11.0/6.0*v3-7.0/6.0*v4+1.0/3.0*v5  # r=0
        return w0 * q0 + w1 * q1 + w2 * q2
