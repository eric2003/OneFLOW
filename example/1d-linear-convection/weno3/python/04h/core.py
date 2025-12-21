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
            
class ComputationalDomain:
    def __init__(self, config, mesh):
        """
        初始化计算域
        :param mesh: Mesh实例（静态网格属性）
        :param config: CfdConfig实例（包含recon_scheme/spatial_order）
        """
        self.config = config
        self.mesh = mesh        
        
        # 核心：根据重建格式动态计算nghosts
        self.nghosts = self._calc_nghosts()
        
        # 基于nghosts推导索引
        self.ist = self.nghosts  # 物理网格起始索引
        self.ied = self.ist + mesh.ncells  # 物理网格结束索引
        self.ntcells = mesh.ncells + 2 * self.nghosts  # 总网格数（含ghost）
       
        print(f"mesh.ncells={mesh.ncells}")
        print(f"self.config.spatial_order={self.config.spatial_order}")
        print(f"self.nghosts={self.nghosts}")
        print(f"self.ist={self.ist}")
        print(f"self.ied={self.ied}")
        
    def _calc_nghosts(self):
        """内部方法：根据重建格式和阶数计算ghost层数量"""
        scheme = self.config.recon_scheme
        order = self.config.spatial_order
        
        if scheme is None:
            raise ValueError("必须先通过 with_reconstruction 设置重建格式！")
        
        # 按格式规则计算nghosts
        if scheme == "eno":
            nghosts = order  # ENO：nghosts = 空间阶数
        elif scheme == "weno":
            nghosts = order // 2 + 1  # WENO：nghosts = 阶数//2 +1
        else:
            raise ValueError(f"未知重建格式 {scheme}，无法计算ghost层！")
        
        # 校验：避免nghosts为0或负数
        if nghosts <= 0:
            raise ValueError(f"计算得到的ghost层数量无效：{nghosts}（阶数{order}，格式{scheme}）")
        return nghosts        
        
    # 可选：提供便捷的索引校验/计算方法，增强复用性
    def is_physical_cell(self, idx):
        """判断索引是否在物理网格范围内"""
        return self.ist <= idx < self.ied
    
    def get_physical_indices(self):
        """返回物理网格的索引范围（可直接用于循环）"""
        return range(self.ist, self.ied)        
        
class Solution:
    def __init__(self, config, domain):
        """
        初始化求解过程中的动态数据
        :param domain: ComputationalDomain实例（用于确定数组尺寸）
        """
        self.domain = domain
        mesh = domain.mesh

        # 界面值和通量（维度依赖mesh.nnodes）
        self.q_face_left = np.zeros(mesh.nnodes)   # 左界面值
        self.q_face_right = np.zeros(mesh.nnodes)  # 右界面值
        self.flux = np.zeros(mesh.nnodes)          # 通量

        # 残差（维度依赖mesh.ncells）
        self.res = np.zeros(mesh.ncells)           # 残差

        # 解数组（维度依赖ntcells，含ghost层）
        self.u = np.zeros(domain.ntcells)          # 当前解
        self.un = np.zeros(domain.ntcells)         # 上一时间步解
        
        self.initialize_from_config(config)

    # 可选：添加数据重置/初始化方法，增强复用性
    def reset_solution(self):
        """重置解数组为初始状态"""
        self.u.fill(0.0)
        self.un.fill(0.0)
        
    def initialize_from_config(self, config):
        """根据配置初始化场"""
        ic = InitialConditionFactory.create(config.ic_type, config)
        ic.apply(self)
        
    def update_old_field(self):
        """更新旧场"""
        #update_oldfield(self.un, self.u)
        self.un[:] = self.u[:]
        
# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, config, mesh):
        self.config = config
        self.domain = ComputationalDomain(config, mesh)
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
