# reconstructor/weno3.py
import numpy as np
from .base import Reconstructor
from cfd_registry import ComponentRegistry, register_component

@register_component('reconstructor', 'weno3')
class Weno3Reconstructor(Reconstructor):
    def compute_face_values(self, q, cfd):
        domain = cfd.domain
        solution = cfd.solution
        self._reconstruct_left_interfaces(domain, q, solution.q_face_left)
        self._reconstruct_right_interfaces(domain, q, solution.q_face_right)
        
    def _reconstruct_left_interfaces(self, domain, u, qL):
        """在每个 i+1/2 界面，计算左单元贡献的 qL (即 u_{i+1/2}^-)"""
        for i in range(domain.ist - 1, domain.ied):
            j = i - (domain.ist - 1)
            v1, v2, v3 = u[i-1], u[i], u[i+1]
            #qL[j] = self._reconstruct_from_right_biased_stencil(v3, v2, v1)
            qL[j] = self._reconstruct_from_left_biased_stencil(v1, v2, v3)
            
    def _reconstruct_right_interfaces(self, domain, u, qR):
        """在每个 i+1/2 界面，计算右单元贡献的 qR (即 u_{i+1/2}^+)"""
        for i in range(domain.ist, domain.ied + 1):
            j = i - domain.ist
            v1, v2, v3 = u[i-1], u[i], u[i+1]
            #qR[j] = self._reconstruct_from_left_biased_stencil(v3, v2, v1)
            qR[j] = self._reconstruct_from_right_biased_stencil(v1, v2, v3)
            
    def _reconstruct_from_left_biased_stencil(self, v1, v2, v3):
        eps = 1e-6
        beta0 = (v2 - v1)**2
        beta1 = (v3 - v2)**2
        d0 = 1.0/3.0
        d1 = 2.0/3.0
        alpha0 = d0 / (eps + beta0)**2
        alpha1 = d1 / (eps + beta1)**2
        alpha = alpha0 + alpha1
        w0 = alpha0 / alpha
        w1 = alpha1 / alpha
        q0 = -1.0/2.0*v1+3.0/2.0*v2  # r=1
        q1 = 1.0/2.0*v2+1.0/2.0*v3  # r=0
        return w0 * q0 + w1 * q1

    def _reconstruct_from_right_biased_stencil(self, v1, v2, v3):
        eps = 1e-6
        beta0 = (v2 - v1)**2
        beta1 = (v3 - v2)**2
        d0 = 2.0/3.0
        d1 = 1.0/3.0
        alpha0 = d0 / (eps + beta0)**2
        alpha1 = d1 / (eps + beta1)**2
        alpha = alpha0 + alpha1
        w0 = alpha0 / alpha
        w1 = alpha1 / alpha
        q0 = 1.0/2.0*v1+1.0/2.0*v2  # r=1
        q1 = 3.0/2.0*v2-1.0/2.0*v3  # r=0
        return w0 * q0 + w1 * q1            

