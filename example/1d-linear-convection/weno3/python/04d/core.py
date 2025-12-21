import numpy as np
import matplotlib.pyplot as plt
import inflect
from abc import ABC, abstractmethod

# 已有 flux 导入
from flux import InviscidFluxCalculator, RusanovFluxCalculator, EngquistOsherFluxCalculator, FluxCalculatorFactory

# 新增 boundary 导入
from boundary import BoundaryCondition, PeriodicBoundary, DirichletBoundary, NeumannBoundary, BoundaryConditionFactory

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

# Initialize reconstruction coefficients for different orders
def init_coef( spatial_order, coef ):
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
    
# ---------------------- 1. 抽象时间推进器基类（统一接口） ----------------------
class TimeIntegrator(ABC):
    """时间推进器抽象基类：定义一维CFD时间推进的核心接口"""
    def __init__(self, cfd):
        self.cfd = cfd  # 持有CFD实例，获取配置/域/求解数据
        self.config = cfd.config
        self.domain = cfd.domain
        self.solution = cfd.solution
        self.residual_calculator = cfd.residual_calculator

    @abstractmethod
    def step(self, dt):
        """
        单次时间步推进（核心接口）
        :param dt: 时间步长
        :return: None
        """
        pass

    # 公共逻辑：复用残差计算、边界条件、数组索引映射
    def compute_residual(self):
        """计算残差（所有RK方法都需要，封装为公共方法）"""
        self.residual_calculator.compute()

    def apply_boundary(self):
        """应用边界条件（公共逻辑）"""
        #boundary(self.solution.u, self.cfd)
        self.cfd.boundary_condition.apply(self.solution.u)
        

    def map_idx(self, i):
        """物理网格索引 → 残差数组索引（公共映射逻辑）"""
        return i - self.domain.ist

# ---------------------- 2. 具体RK时间推进器实现（复用公共逻辑） ----------------------
class RK1Integrator(TimeIntegrator):
    """1阶显式欧拉（RK1）"""
    def step(self, dt):
        # RK1核心逻辑：u = u + dt * res
        self.compute_residual()  # 复用公共残差计算
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] += dt * self.solution.res[j]
        self.apply_boundary()  # 复用公共边界条件
        self.solution.update_old_field()  # 同步old field

class RK2Integrator(TimeIntegrator):
    """2阶Heun方法（RK2）"""
    def step(self, dt):
        # 阶段1：预测步
        self.compute_residual()
        u_pred = self.solution.u.copy()  # 保存预测值
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u_pred[i] += dt * self.solution.res[j]
        self.solution.u[:] = u_pred  # 更新预测值
        self.apply_boundary()

        # 阶段2：校正步
        self.compute_residual()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] = 0.5 * self.solution.un[i] + 0.5 * self.solution.u[i] + 0.5 * dt * self.solution.res[j]
        self.apply_boundary()
        self.solution.update_old_field()

class RK3Integrator(TimeIntegrator):
    """3阶SSPRK3（强稳定保号RK3）"""
    def step(self, dt):
        # 阶段1
        self.compute_residual()
        u1 = self.solution.u.copy()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u1[i] += dt * self.solution.res[j]
        self.solution.u[:] = u1
        self.apply_boundary()

        # 阶段2
        self.compute_residual()
        u2 = self.solution.u.copy()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u2[i] = 0.75 * self.solution.un[i] + 0.25 * self.solution.u[i] + 0.25 * dt * self.solution.res[j]
        self.solution.u[:] = u2
        self.apply_boundary()

        # 阶段3
        self.compute_residual()
        c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] = c1 * self.solution.un[i] + c2 * self.solution.u[i] + c3 * dt * self.solution.res[j]
        self.apply_boundary()
        self.solution.update_old_field()

# ---------------------- 3. 时间推进器工厂（统一创建逻辑） ----------------------
class TimeIntegratorFactory:
    """时间推进器工厂：根据配置创建对应RK实例"""
    @staticmethod
    def create(cfd):
        rk_order = cfd.config.rk_order
        integrator_mapping = {
            1: RK1Integrator,
            2: RK2Integrator,
            3: RK3Integrator,
            # 新增RK4只需：4: RK4Integrator
        }
        if rk_order not in integrator_mapping:
            raise ValueError(f"不支持的RK阶数：{rk_order}（可选：{list(integrator_mapping.keys())}）")
        return integrator_mapping[rk_order](cfd)    

# Mesh class: defines computational grid
class Mesh:
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 2.0
        self.ncells = 40
        self.nnodes = self.ncells + 1
        self.nx = self.ncells
        self.x = np.zeros(self.nnodes)
        self.xcc = np.zeros(self.ncells)
        self.init_mesh()
        
    def init_mesh(self):
        self.L = self.xmax - self.xmin
        self.dx = self.L / self.ncells
        
        # Generate node coordinates
        for i in range(self.nnodes):
            self.x[i] = self.xmin + i * self.dx
            
        # Generate cell center coordinates
        for i in range(self.ncells):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i+1])

class Reconstructor(ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def reconstruct(self, q, cfd):
        pass

class EnoReconstructor(Reconstructor):
    def __init__(self, spatial_order, ntcells):
        self.spatial_order = spatial_order
        self.ntcells = ntcells
        # Stencil selection arrays
        self.lmc = np.zeros(self.ntcells, dtype=int)
        
        # Reconstruction coefficients and divided differences
        self.coef = np.zeros((spatial_order + 1, spatial_order))
        self.dd = np.zeros((spatial_order, self.ntcells))
        
        init_coef(self.spatial_order, self.coef)

    def reconstruct(self, q, cfd):
        """ENO reconstruction of interface values"""
        
        # Choose stencil by ENO method based on smoothest polynomial
        self.dd[0, :] = q
        
        # Compute divided differences
        for m in range(1, self.spatial_order):
            for j in range(self.ntcells-m):
                self.dd[m, j] = self.dd[m-1, j+1] - self.dd[m-1, j]
                
        domain = cfd.domain
        solution = cfd.solution
                
        # Select left-biased stencil for each node
        for i in range(domain.ist-1,domain.ied+1):
            self.lmc[i] = i
            for m in range(1, self.spatial_order):
                if abs(self.dd[m, self.lmc[i]-1]) < abs(self.dd[m, self.lmc[i]]):
                    self.lmc[i] -= 1
                    
        # Reconstruct values at cell interfaces (j+1/2)
        for i in range(domain.ist,domain.ied+1):
            j = i - domain.ist
            k1 = self.lmc[i-1]
            k2 = self.lmc[i  ]
            r1 = i-1 - k1
            r2 = i   - k2
            #print(f"i,k1,k2,r1,r2={i,k1,k2,r1,r2}")
            solution.q_face_left[j] = 0
            solution.q_face_right[j] = 0
            for m in range(self.spatial_order):
                solution.q_face_left[j] += q[k1 + m] * self.coef[r1+1, m]
                solution.q_face_right[j] += q[k2 + m] * self.coef[r2, m]

    
class WenoReconstructor(Reconstructor):
    def reconstruct(self, q, cfd):
        # WENO (Weighted Essentially Non-Oscillatory) reconstruction
        # Reconstruct values at cell interfaces (j+1/2)
        domain = cfd.domain
        solution = cfd.solution
        self.weno3L( domain, q, solution.q_face_left )
        self.weno3R( domain, q, solution.q_face_right )
        
    # 3rd-order WENO reconstruction for left interface with periodic boundary
    def weno3L(self, domain, u, f):
        # i: ist-1, ist, ..., ied-1
        # j: 0, 1, ..., nx
        for i in range(domain.ist - 1, domain.ied):
            j = i - ( domain.ist - 1 )
            v1 = u[i-1]
            v2 = u[i  ]
            v3 = u[i+1]
            #f[j] = self.wc3L(v1,v2,v3)
            f[j] = self.wc3R(v3,v2,v1)
            
    # 3rd-order WENO reconstruction for right interface with periodic boundary
    def weno3R(self, domain, u, f):
        # i: ist, ist+1, ..., ied, ied+1
        # j: 0, 1, ..., nx
        for i in range(domain.ist, domain.ied + 1):
            j = i - domain.ist
            v1 = u[i-1]
            v2 = u[i  ]
            v3 = u[i+1]
            #f[j] = self.wc3R(v1,v2,v3)
            f[j] = self.wc3L(v3,v2,v1)
            
    def wc3L(self,v1,v2,v3):
        """WENO-3 nonlinear weights for left-biased stencil"""
        eps = 1.0e-6

        # Smoothness indicators
        s0 = (v3-v2)**2
        s1 = (v2-v1)**2

        # Compute nonlinear weights w0, w1
        d0 = 2.0/3.0
        d1 = 1.0/3.0
        
        c0 = d0 / ( (eps+s0)**2 )
        c1 = d1 / ( (eps+s1)**2 )

        w0 = c0 / ( c0 + c1 )
        w1 = c1 / ( c0 + c1 )

        # Candidate stencils
        q0 =  0.5 * v2 + 0.5 * v3
        q1 = -0.5 * v1 + 1.5 * v2

        # Reconstructed value at interface
        f = ( w0*q0 + w1*q1 )

        return f

    def wc3R(self,v1,v2,v3):
        """WENO-3 nonlinear weights for right-biased stencil"""
        eps = 1.0e-6

        # Smoothness indicators
        s0 = (v3-v2)**2
        s1 = (v2-v1)**2

        # Compute nonlinear weights w0, w1
        d0 = 1.0/3.0
        d1 = 2.0/3.0
        
        c0 = d0 / ( (eps+s0)**2 )
        c1 = d1 / ( (eps+s1)**2 )

        w0 = c0 / ( c0 + c1 )
        w1 = c1 / ( c0 + c1 )

        # Candidate stencils
        q0 = 1.5 * v2 - 0.5 * v3
        q1 = 0.5 * v1 + 0.5 * v2

        # Reconstructed value at interface
        f = ( w0*q0 + w1*q1 )

        return f            
        
class ReconstructorFactory:
    """重建器工厂：根据配置自动创建对应类型的重建器实例"""
    @staticmethod
    def create(config, domain):
        """
        静态工厂方法（无需实例化工厂类）
        :param config: CfdConfig实例（含recon_scheme/spatial_order）
        :param domain: ComputationalDomain实例（含ntcells）
        :return: Reconstructor子类实例
        """
        scheme = config.recon_scheme.lower()
        if scheme == "eno":
            # ENO需要空间阶数和总网格数
            return EnoReconstructor(
                spatial_order=config.spatial_order,
                ntcells=domain.ntcells
            )
        elif scheme == "weno":
            # WENO无需额外参数（可根据需求扩展，如传入order）
            return WenoReconstructor()
        else:
            raise ValueError(f"不支持的重建格式：{scheme}（仅支持 eno/weno）")        
        
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
        
class InitialCondition(ABC):
    """初始条件基类"""
    def __init__(self, config):
        self.config = config

    @abstractmethod
    def apply(self, solution):
        """将初始条件应用到 solution 的内部区域"""
        pass

    @abstractmethod
    def evaluate_at(self, x):
        """纯数学函数：给定 x，返回 u(x)，不涉及网格或边界"""
        pass

    def _apply_to_interior(self, solution, values):
        domain = solution.domain
        for i in range(domain.ist, domain.ied):
            j = i - domain.ist
            solution.u[i] = values[j]


class StepFunctionIC(InitialCondition):
    def evaluate_at(self, x):
        u0 = np.ones_like(x)
        mask = (x >= 0.5) & (x <= 1.0)
        u0[mask] = 2.0
        return u0

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)


class SineWaveIC(InitialCondition):
    def evaluate_at(self, x):
        L = self.config.get("domain_length", 2.0)
        return np.sin(2 * np.pi * x / L)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)


class GaussianPulseIC(InitialCondition):
    def evaluate_at(self, x):
        center = self.config.get("pulse_center", 0.5)
        width = self.config.get("pulse_width", 0.1)
        return np.exp(-((x - center) / width) ** 2)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)


class InitialConditionFactory:
    _registry = {
        'step': StepFunctionIC,
        'sin': SineWaveIC,
        'gaussian': GaussianPulseIC,
    }

    @classmethod
    def create(cls, ic_type, config):
        if ic_type not in cls._registry:
            raise ValueError(f"未知的初始条件类型: {ic_type}（支持: {list(cls._registry.keys())}）")
        return cls._registry[ic_type](config)

    @classmethod
    def register(cls, name, ic_class):
        cls._registry[name] = ic_class


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
