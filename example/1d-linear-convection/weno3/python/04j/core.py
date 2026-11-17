import numpy as np
import matplotlib.pyplot as plt
import inflect
from abc import ABC, abstractmethod

# Flux
from flux import InviscidFluxCalculator, RusanovFluxCalculator, EngquistOsherFluxCalculator, FluxCalculatorFactory

# Boundary
from boundary import BoundaryCondition, PeriodicBoundary, DirichletBoundary, NeumannBoundary, BoundaryConditionFactory

# Time integration
from time_integration import TimeIntegrator, RK1Integrator, RK2Integrator, RK3Integrator, TimeIntegratorFactory

# Mesh 👈 新增这一行
from mesh import Mesh

from reconstructor import Reconstructor, EnoReconstructor, WenoReconstructor, ReconstructorFactory, init_coef

from initial_condition import InitialCondition, StepFunctionIC, SineWaveIC, GaussianPulseIC, InitialConditionFactory

from domain import Domain
from solution import Solution  # 👈 新增

# ---------------------- 4. 残差计算器（封装完整残差计算逻辑） ----------------------
class ResidualCalculator:
    """残差计算器：封装「重建→通量→散度」完整流程"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.config = cfd.config
        self.domain = cfd.domain
        self.solution = cfd.solution
        self.mesh = self.domain.mesh
        self.reconstructor = self.cfd.reconstructor
        
        # 初始化通量计算器（工厂创建）
        self.flux_calculator = FluxCalculatorFactory.create(cfd)

    def compute(self):
        """计算完整残差（对外唯一接口）"""
        # 步骤1：界面重建（调用外部重建函数，保持兼容）
        self._reconstruct()
        
        # 步骤2：计算无粘通量
        self._compute_inviscid_flux()
        
        # 步骤3：计算通量散度（残差核心）
        self._compute_flux_divergence()

    def _reconstruct(self):
        """私有方法：界面值重建"""
        self.reconstructor.reconstruct(self.solution.u, self.cfd)

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
        # 向量化计算：残差[i] = -(flux[i+1] - flux[i])/dx
        for i in range(self.mesh.ncells):
            solution.res[i] = -(solution.flux[i+1] - solution.flux[i]) / self.mesh.dx
            
class CfdConfig:
    def __init__(self):
        self.ic_type = "step"
        self.recon_scheme = "eno"  # 0=ENO, 1=WENO
        self.flux_type = 0     # 0=Rusanov, 1=Engquist-Osher
        self.rk_order = 1
        self.wave_speed = 1.0
        self.final_time = 0.625
        self.dt = 0.025
        
        self.boundary_type = "periodic"
        self.left_boundary_value = 1.0  # Dirichlet左边界值
        self.right_boundary_value = 2.0 # Dirichlet右边界值
        
    def with_reconstruction(self, scheme, order=None):
        """专用配置：重建方案（链式调用）"""
        self.recon_scheme = scheme.lower()  # 统一小写，避免大小写问题
        
        # 智能默认阶数
        if order is not None:
            self.spatial_order = order
        else:
            if self.recon_scheme == "weno":
                self.spatial_order = 5
            elif self.recon_scheme == "eno":
                self.spatial_order = 3  # ENO默认3阶
            else:
                raise ValueError(f"不支持的重建格式：{scheme}（仅支持 eno/weno）")
        
        return self
        
    def with_boundary(self, bc_type, left_value=None, right_value=None):
        self.boundary_type = bc_type
        if left_value is not None:
            self.left_boundary_value = left_value
        if right_value is not None:
            self.right_boundary_value = right_value
        return self        
            
# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, config, mesh):
        self.config = config
        self.domain = Domain(config, mesh)
        self.solution = Solution(config, self.domain)
        self.reconstructor = ReconstructorFactory.create(config, self.domain)
        self.residual_calculator = ResidualCalculator(self)        
        self.integrator = TimeIntegratorFactory.create(self)
        self.boundary_condition = BoundaryConditionFactory.create(self)
    
    def exact_solution(self):
        """通用对流问题的解析解：u(x, T) = u0(x - c*T)，周期边界"""
        x = self.domain.mesh.xcc
        T = self.config.final_time
        c = self.config.wave_speed
        L = self.domain.mesh.L
 
        # 周期平移：确保在 [x0, x0 + L) 内
        x_shifted = (x - c * T + L) % L

        # 获取 IC 实例并评估
        ic = InitialConditionFactory.create(self.config.ic_type, self.config)
        return ic.evaluate_at(x_shifted)
    
    def run(self):
        # 应用初始边界条件并同步 old field
        self.boundary_condition.apply(self.solution.u)
        self.solution.update_old_field()       
        
        t = 0.0
        dt_old = self.config.dt
        dt = dt_old
        
        domain = self.domain
        solution = self.solution
        config = self.config
        
        while t < config.final_time:
            if t + dt > config.final_time:
                dt = config.final_time - t
            config.dt = dt  # temporary adjustment for last step
            #runge_kutta(self)
            self.integrator.step(dt)
            t += dt
        config.dt = dt_old
        
        # 整理标准化结果
        u_numerical = self.solution.u[self.domain.ist:self.domain.ied].copy()
        self.result = {
            "x": domain.mesh.xcc,
            "numerical": u_numerical,
            "analytical": self.exact_solution(),
            "config": {
                "scheme": self.config.recon_scheme,
                "order": self.config.spatial_order,
                "rk_order": self.config.rk_order,
                "final_time": self.config.final_time
            }
        }

        return u_numerical
