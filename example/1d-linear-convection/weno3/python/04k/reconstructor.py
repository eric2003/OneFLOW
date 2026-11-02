# reconstructor.py
import numpy as np
from abc import ABC, abstractmethod

# ---------------------- 1. 重构系数初始化函数 ----------------------
def init_coef(spatial_order, coef):
    """Initialize reconstruction coefficients for different spatial orders."""
    if spatial_order == 1:
        coef[0] = [1.0]
        coef[1] = [1.0]
    elif spatial_order == 2:
        coef[0] = [3.0/2.0, -1.0/2.0]
        coef[1] = [1.0/2.0,  1.0/2.0]
        coef[2] = [-1.0/2.0, 3.0/2.0]
    elif spatial_order == 3:
        coef[0] = [ 11.0/6.0, -7.0/6.0,  1.0/3.0 ]
        coef[1] = [  1.0/3.0,  5.0/6.0, -1.0/6.0 ]
        coef[2] = [ -1.0/6.0,  5.0/6.0,  1.0/3.0 ]
        coef[3] = [  1.0/3.0, -7.0/6.0, 11.0/6.0 ]
    elif spatial_order == 4:
        coef[0] = [ 25.0/12.0, -23.0/12.0,  13.0/12.0,  -1.0/4.0 ]
        coef[1] = [   1.0/4.0,  13.0/12.0,  -5.0/12.0,  1.0/12.0 ]
        coef[2] = [ -1.0/12.0,   7.0/12.0,   7.0/12.0, -1.0/12.0 ]
        coef[3] = [  1.0/12.0,  -5.0/12.0,  13.0/12.0,   1.0/4.0 ]
        coef[4] = [  -1.0/4.0,  13.0/12.0, -23.0/12.0, 25.0/12.0 ]
    elif spatial_order == 5:
        coef[0] = [ 137.0/60.0, -163.0/60.0, 137.0/60.0,  -21.0/20.0,    1.0/5.0 ]
        coef[1] = [    1.0/5.0,   77.0/60.0, -43.0/60.0,   17.0/60.0,  -1.0/20.0 ]
        coef[2] = [  -1.0/20.0,    9.0/20.0,  47.0/60.0,  -13.0/60.0,   1.0/30.0 ]
        coef[3] = [   1.0/30.0,  -13.0/60.0,  47.0/60.0,    9.0/20.0,  -1.0/20.0 ]
        coef[4] = [  -1.0/20.0,   17.0/60.0, -43.0/60.0,   77.0/60.0,    1.0/5.0 ]
        coef[5] = [    1.0/5.0,  -21.0/20.0, 137.0/60.0, -163.0/60.0, 137.0/60.0 ]
    elif spatial_order == 6:
        coef[0] = [ 49.0/20.0, -71.0/20.0,   79.0/20.0, -163.0/60.0,  31.0/30.0,  -1.0/6.0 ]
        coef[1] = [   1.0/6.0,  29.0/20.0,  -21.0/20.0,   37.0/60.0, -13.0/60.0,  1.0/30.0 ]
        coef[2] = [ -1.0/30.0,  11.0/30.0,   19.0/20.0,  -23.0/60.0,   7.0/60.0, -1.0/60.0 ]
        coef[3] = [  1.0/60.0,  -2.0/15.0,   37.0/60.0,   37.0/60.0,  -2.0/15.0,  1.0/60.0 ]
        coef[4] = [ -1.0/60.0,   7.0/60.0,  -23.0/60.0,   19.0/20.0,  11.0/30.0, -1.0/30.0 ]
        coef[5] = [  1.0/30.0, -13.0/60.0,   37.0/60.0,  -21.0/20.0,  29.0/20.0,   1.0/6.0 ]
        coef[6] = [  -1.0/6.0,  31.0/30.0, -163.0/60.0,   79.0/20.0, -71.0/20.0, 49.0/20.0 ]
    elif spatial_order == 7:
        coef[0] = [ 363.0/140.0, -617.0/140.0,  853.0/140.0, -2341.0/420.0,  667.0/210.0,   -43.0/42.0,     1.0/7.0 ]
        coef[1] = [     1.0/7.0,  223.0/140.0, -197.0/140.0,   153.0/140.0, -241.0/420.0,   37.0/210.0,   -1.0/42.0 ]
        coef[2] = [   -1.0/42.0,    13.0/42.0,  153.0/140.0,  -241.0/420.0,  109.0/420.0,  -31.0/420.0,   1.0/105.0 ]
        coef[3] = [   1.0/105.0,  -19.0/210.0,  107.0/210.0,   319.0/420.0, -101.0/420.0,     5.0/84.0,  -1.0/140.0 ]
        coef[4] = [  -1.0/140.0,     5.0/84.0, -101.0/420.0,   319.0/420.0,  107.0/210.0,  -19.0/210.0,   1.0/105.0 ]
        coef[5] = [   1.0/105.0,  -31.0/420.0,  109.0/420.0,  -241.0/420.0,  153.0/140.0,    13.0/42.0,   -1.0/42.0 ]
        coef[6] = [   -1.0/42.0,   37.0/210.0, -241.0/420.0,   153.0/140.0, -197.0/140.0,  223.0/140.0,     1.0/7.0 ]
        coef[7] = [     1.0/7.0,   -43.0/42.0,  667.0/210.0, -2341.0/420.0,  853.0/140.0, -617.0/140.0, 363.0/140.0 ]


# ---------------------- 2. 抽象重构器基类 ----------------------
class Reconstructor(ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def reconstruct(self, q, cfd):
        pass


# ---------------------- 3. ENO 重构器 ----------------------
class EnoReconstructor(Reconstructor):
    def __init__(self, spatial_order, ntcells):
        self.spatial_order = spatial_order
        self.ntcells = ntcells
        self.lmc = np.zeros(self.ntcells, dtype=int)
        self.coef = np.zeros((spatial_order + 1, spatial_order))
        self.dd = np.zeros((spatial_order, self.ntcells))
        init_coef(self.spatial_order, self.coef)

    def reconstruct(self, q, cfd):
        """ENO reconstruction of interface values"""
        self.dd[0, :] = q
        for m in range(1, self.spatial_order):
            for j in range(self.ntcells - m):
                self.dd[m, j] = self.dd[m-1, j+1] - self.dd[m-1, j]
                
        domain = cfd.domain
        solution = cfd.solution
                
        for i in range(domain.ist - 1, domain.ied + 1):
            self.lmc[i] = i
            for m in range(1, self.spatial_order):
                if abs(self.dd[m, self.lmc[i] - 1]) < abs(self.dd[m, self.lmc[i]]):
                    self.lmc[i] -= 1
                    
        for i in range(domain.ist, domain.ied + 1):
            j = i - domain.ist
            k1 = self.lmc[i - 1]
            k2 = self.lmc[i]
            r1 = i - 1 - k1
            r2 = i - k2
            solution.q_face_left[j] = 0.0
            solution.q_face_right[j] = 0.0
            for m in range(self.spatial_order):
                solution.q_face_left[j] += q[k1 + m] * self.coef[r1 + 1, m]
                solution.q_face_right[j] += q[k2 + m] * self.coef[r2, m]


# ---------------------- 4. WENO 重构器（3阶） ----------------------
class WenoReconstructor(Reconstructor):
    def reconstruct(self, q, cfd):
        domain = cfd.domain
        solution = cfd.solution
        self.weno3L(domain, q, solution.q_face_left)
        self.weno3R(domain, q, solution.q_face_right)
        
    def weno3L(self, domain, u, f):
        for i in range(domain.ist - 1, domain.ied):
            j = i - (domain.ist - 1)
            v1 = u[i - 1]
            v2 = u[i]
            v3 = u[i + 1]
            f[j] = self.wc3R(v3, v2, v1)
            
    def weno3R(self, domain, u, f):
        for i in range(domain.ist, domain.ied + 1):
            j = i - domain.ist
            v1 = u[i - 1]
            v2 = u[i]
            v3 = u[i + 1]
            f[j] = self.wc3L(v3, v2, v1)
            
    def wc3L(self, v1, v2, v3):
        eps = 1.0e-6
        s0 = (v3 - v2)**2
        s1 = (v2 - v1)**2
        d0, d1 = 2.0/3.0, 1.0/3.0
        c0 = d0 / ((eps + s0)**2)
        c1 = d1 / ((eps + s1)**2)
        w0 = c0 / (c0 + c1)
        w1 = c1 / (c0 + c1)
        q0 = 0.5 * v2 + 0.5 * v3
        q1 = -0.5 * v1 + 1.5 * v2
        return w0 * q0 + w1 * q1

    def wc3R(self, v1, v2, v3):
        eps = 1.0e-6
        s0 = (v3 - v2)**2
        s1 = (v2 - v1)**2
        d0, d1 = 1.0/3.0, 2.0/3.0
        c0 = d0 / ((eps + s0)**2)
        c1 = d1 / ((eps + s1)**2)
        w0 = c0 / (c0 + c1)
        w1 = c1 / (c0 + c1)
        q0 = 1.5 * v2 - 0.5 * v3
        q1 = 0.5 * v1 + 0.5 * v2
        return w0 * q0 + w1 * q1


# ---------------------- 5. 重构器工厂 ----------------------
class ReconstructorFactory:
    """重建器工厂：根据配置自动创建对应类型的重建器实例"""
    @staticmethod
    def create(config, domain):
        scheme = config.recon_scheme.lower()
        if scheme == "eno":
            return EnoReconstructor(
                spatial_order=config.spatial_order,
                ntcells=domain.ntcells
            )
        elif scheme == "weno":
            return WenoReconstructor()
        else:
            raise ValueError(f"不支持的重建格式：{scheme}（仅支持 eno/weno）")