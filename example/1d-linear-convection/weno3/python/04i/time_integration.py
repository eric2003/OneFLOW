# time_integration.py
from abc import ABC, abstractmethod

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
        self.cfd.boundary_condition.apply(self.solution.u)

    def map_idx(self, i):
        """物理网格索引 → 残差数组索引（公共映射逻辑）"""
        return i - self.domain.ist

# ---------------------- 2. 具体RK时间推进器实现（复用公共逻辑） ----------------------
class RK1Integrator(TimeIntegrator):
    """1阶显式欧拉（RK1）"""
    def step(self, dt):
        self.compute_residual()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] += dt * self.solution.res[j]
        self.apply_boundary()
        self.solution.update_old_field()

class RK2Integrator(TimeIntegrator):
    """2阶Heun方法（RK2）"""
    def step(self, dt):
        # 阶段1：预测步
        self.compute_residual()
        u_pred = self.solution.u.copy()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            u_pred[i] += dt * self.solution.res[j]
        self.solution.u[:] = u_pred
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