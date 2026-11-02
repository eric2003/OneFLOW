import numpy as np
import matplotlib.pyplot as plt
import inflect
from abc import ABC, abstractmethod

def initial_condition(x):
    """Initial condition: step function from 1.0 to 2.0 in [0.5, 1.0]"""
    u0 = np.zeros_like(x)
    for i in range(len(x)):
        if 0.5 <= x[i] <= 1.0:
            u0[i] = 2.0
        else:
            u0[i] = 1.0
    return u0

def analytical_solution(x, t, a, L):
    """Analytical solution with periodic boundary conditions"""
    x_shifted = (x - a * t + L) % L
    return initial_condition(x_shifted)

# Compute residual (flux divergence) for all cells
def residual(q, cfd):
    reconstruction(q, cfd)
    solution = cfd.solution  
    mesh = cfd.domain.mesh    
    inviscid_flux(solution.q_face_left, solution.q_face_right, solution.flux, cfd)
    for i in range(mesh.ncells):
        solution.res[i] = -(solution.flux[i+1] - solution.flux[i]) / mesh.dx
        
# Choose reconstruction method based on solver setting
def reconstruction(q, cfd):
    cfd.reconstructor.reconstruct(q, cfd)

# Compute inviscid flux using selected Riemann solver
def inviscid_flux(q_face_left, q_face_right, flux, cfd):
    if cfd.config.flux_type == 0:
        rusanov_flux(q_face_left, q_face_right, flux, cfd)
    else:
        engquist_osher_flux(q_face_left, q_face_right, flux, cfd)
        
# --------------------------------------------------------------------------- #
# Numerical fluxes
# --------------------------------------------------------------------------- #
def rusanov_flux(q_face_left, q_face_right, flux, cfd):
    """Rusanov (local Lax-Friedrichs) flux"""
    mesh = cfd.domain.mesh

    for i in range(mesh.nnodes):
        u_L = q_face_left[i]
        u_R = q_face_right[i]
        c_L = cfd.config.wave_speed
        c_R = cfd.config.wave_speed
        F_L = c_L * u_L  # Flux from left state
        F_R = c_R * u_R  # Flux from right state
        Smax = max(abs(c_L),abs(c_R)) # Maximum wave speed
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L)
        
def engquist_osher_flux(q_face_left, q_face_right, flux, cfd):
    """Engquist-Osher flux for linear convection"""
    for i in range(cfd.nnodes):
        c = cfd.config.wave_speed
        cp = 0.5*(c + abs(c))
        cm = 0.5*(c - abs(c))
        u_L = q_face_left[i]
        u_R = q_face_right[i]
       
        flux[i] = cp * u_L + cm * u_R

def periodic_boundary(u, cfd):
    """Apply periodic boundary conditions"""
    domain = cfd.domain
    # Left ghost cells = right interior cells
    for ig in range(domain.nghosts):
        u[domain.ist - 1 - ig] = u[domain.ied - 1 - ig]
    
    # Right ghost cells = left interior cells
    for ig in range(domain.nghosts):
        u[domain.ied + ig] = u[domain.ist + ig]
    
    
# Apply periodic boundary conditions
def boundary(u, cfd):
    periodic_boundary(u, cfd)

# Copy current solution to old solution array
def update_oldfield(qn, q):
    qn[:] = q[:]

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

# Select Runge-Kutta time integration scheme
def runge_kutta(cfd):
    rk_order = cfd.config.rk_order
    if rk_order == 1:
        runge_kutta_1(cfd)
    elif rk_order == 2:
        runge_kutta_2(cfd)
    else:
        runge_kutta_3(cfd)
    
# 1st-order explicit Euler time integration
def runge_kutta_1(cfd):
    dt = cfd.config.dt
    domain = cfd.domain
    solution = cfd.solution
    
    residual(solution.u, cfd)
    for i in range(domain.ist,domain.ied):
        j = i - domain.ist
        solution.u[i] = solution.u[i] + dt * solution.res[j]
    boundary(solution.u, cfd)
    update_oldfield(solution.un, solution.u)
    
# 2nd-order Runge-Kutta (Heun's method) time integration
def runge_kutta_2(cfd):
    dt = cfd.config.dt
    domain = cfd.domain
    solution = cfd.solution
    
    residual(solution.u, cfd)
    
    for i in range(domain.ist,domain.ied):
        j = i - domain.ist
        solution.u[i] = solution.u[i] + dt * solution.res[j]
    boundary(solution.u, cfd)
    
    residual(solution.u, cfd)
   
    for i in range(domain.ist,domain.ied):
        j = i - domain.ist
        solution.u[i] = 0.5 * solution.un[i] + 0.5 * solution.u[i] + 0.5 * dt * solution.res[j]
    boundary(solution.u, cfd)
    update_oldfield(solution.un, solution.u)

# 3rd-order Runge-Kutta (SSPRK3) time integration
def runge_kutta_3(cfd):
    dt = cfd.config.dt
    domain = cfd.domain
    solution = cfd.solution
    
    residual(solution.u, cfd)
    for i in range(domain.ist,domain.ied):
        j = i - domain.ist
        solution.u[i] = solution.u[i] + dt * solution.res[j]
    boundary(solution.u, cfd)
    
    residual(solution.u, cfd)
    for i in range(domain.ist,domain.ied):
        j = i - domain.ist
        solution.u[i] = 0.75 * solution.un[i] + 0.25 * solution.u[i] + 0.25 * dt * solution.res[j]
    boundary(solution.u, cfd)
    
    residual(solution.u, cfd)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in range(domain.ist,domain.ied):
        j = i - domain.ist
        solution.u[i] = c1 * solution.un[i] + c2 * solution.u[i] + c3 * dt * solution.res[j]
    boundary(solution.u, cfd)
    update_oldfield(solution.un, solution.u)    

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
        #print(f"EnoReconstructor:reconstruct")
        
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
        self.recon_scheme = "eno"  # 0=ENO, 1=WENO
        self.flux_type = 0     # 0=Rusanov, 1=Engquist-Osher
        self.rk_order = 1
        self.wave_speed = 1.0
        self.final_time = 0.625
        self.dt = 0.025
        
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
    def __init__(self, domain):
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

    # 可选：添加数据重置/初始化方法，增强复用性
    def reset_solution(self):
        """重置解数组为初始状态"""
        self.u.fill(0.0)
        self.un.fill(0.0)
        
# Cfd class: main data structure containing all CFD data
class Cfd:
    def __init__(self, config, mesh):
        self.config = config
        domain = ComputationalDomain(config, mesh)
        self.domain = domain
        self.solution = Solution(domain)  # 核心求解数据
        self.reconstructor = ReconstructorFactory.create(config, domain)
        
    def init_field(self):
        domain = self.domain
        solution = self.solution
        for i in range(domain.ist, domain.ied):
            j = i - domain.ist
            if 0.5 <= domain.mesh.xcc[j] <= 1.0:
                solution.u[i] = 2.0
            else:
                solution.u[i] = 1.0
        boundary(solution.u, self)
        update_oldfield(solution.un, solution.u)
        
    def initial_condition(self, x):
        """Initial condition: step function from 1.0 to 2.0 in [0.5, 1.0]"""
        u0 = np.zeros_like(x)
        for i in range(len(x)):
            if 0.5 <= x[i] <= 1.0:
                u0[i] = 2.0
            else:
                u0[i] = 1.0
        return u0
        
    def init_field_new(self):
        domain = self.domain
        solution = self.solution
        sol = self.initial_condition(domain.mesh.xcc)
        for i in range(domain.ist, domain.ied):
            j = i - domain.ist
            solution.u[i] = sol[j]
        boundary(solution.u, self)
        update_oldfield(solution.un, solution.u)
        
    def exact_solution(self):
        """Analytical solution with periodic boundary conditions"""
        x = self.mesh.xcc
        T = self.config.final_time
        c = self.config.wave_speed
        L = self.mesh.L
        x_shifted = (x - c * T + L) % L
        return initial_condition(x_shifted)
    
    def run(self):
        #self.init_field()
        self.init_field_new()
        
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
            runge_kutta(self)
            t += dt
        config.dt = dt_old
        return solution.u[domain.ist:domain.ied].copy()


# Perform ENO-WENO comparative analysis
def performEnoWenoAnalysis():
    mesh = Mesh()
    
    config_eno3 = (CfdConfig()
                    .with_reconstruction("eno",3))
                
    config_eno3.rk_order = 1
    config_eno3.dt = 0.0025

    cfd_eno3 = Cfd(config_eno3, mesh)
    
    
    u_list = []
    # ENO
    u_eno = cfd_eno3.run()
    u_list.append(u_eno)
    
    # WENO
    config_weno3 = (CfdConfig()
                    .with_reconstruction("weno",3))
                
    config_weno3.rk_order = 1
    config_weno3.dt = 0.0025

    cfd_weno3 = Cfd(config_weno3, mesh)
    u_weno = cfd_weno3.run()
    u_list.append(u_weno)
    
    u_analytical = analytical_solution(mesh.xcc, config_weno3.final_time, config_weno3.wave_speed, mesh.L)
    plot_EnoWeno_Analysis(config_weno3, mesh.xcc, u_list, u_analytical)
    
# Plot ENO-WENO comparison results
def plot_EnoWeno_Analysis(config, xcc, u_list, u_analytical):
    # Define line styles with different colors and markers
    styles = [
        {'color': 'black', 'linestyle': '-', 'marker': 'o'},
        {'color': 'blue', 'linestyle': '--', 'marker': 's'},
        {'color': 'black', 'linestyle': '-', 'marker': '^'},
        {'color': 'blue', 'linestyle': '--', 'marker': 'v'},
        {'color': 'black', 'linestyle': '-', 'marker': '<'},
        {'color': 'blue', 'linestyle': '--', 'marker': '>'},
        {'color': 'black', 'linestyle': '-', 'marker': 'D'},
    ]

    n = len(u_list)
    num_styles = len(styles)

    p = inflect.engine()
    rk_str = p.ordinal(config.rk_order)    
    
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
    plt.title(f'1D Convection Equation at t = {config.final_time:.3f} using 3rd-order ENO&WENO and {rk_str}-order Runge-Kutta methods')
    for i in range(0, n):
        if i == 0:
            lable = 'Numerical (Rusanov)ENO3'
        else:
            lable = 'Numerical (Rusanov)WENO3'
        style = styles[i % num_styles]
        plt.plot(xcc, u_list[i], marker=style['marker'], markerfacecolor='none', linestyle=style['linestyle'], color=style['color'], \
                markersize=5, linewidth=0.5, alpha=1.0, label=f'{lable}')
    plt.plot(xcc, u_analytical, 'r--', label='Analytical')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    plt.show()    
    
if __name__ == "__main__":
    performEnoWenoAnalysis()