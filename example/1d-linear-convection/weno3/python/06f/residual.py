# residual.py

from flux.factory import FluxCalculatorFactory
from reconstructor import ReconstructorFactory

class ResidualCalculator:
    """残差计算器：封装「重建→通量→散度」完整流程"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.config = cfd.config
        self.domain = cfd.domain
        self.solution = cfd.solution
        self.mesh = self.domain.mesh
        
        self.reconstructor = ReconstructorFactory.create(self.config, self.domain)
        self.flux_calculator = FluxCalculatorFactory.create(cfd)
        
    def compute(self):
        """计算完整残差（对外唯一接口）"""
        self._compute_face_values()
        self._compute_inviscid_flux()
        self._compute_flux_divergence()

    def _compute_face_values(self):
        """私有方法：界面值重建"""
        self.reconstructor.compute_face_values(self.solution.u, self.cfd)

    def _compute_inviscid_flux(self):
        """私有方法：计算无粘通量"""
        self.flux_calculator.compute(
            self.solution.q_face_left,
            self.solution.q_face_right,
            self.solution.flux
        )

    def _compute_flux_divergence(self):
        """私有方法：计算通量散度（残差 = -dF/dx）"""
        solution = self.solution
        for i in range(self.mesh.ncells):
            solution.res[i] = -(solution.flux[i+1] - solution.flux[i]) / self.mesh.dx