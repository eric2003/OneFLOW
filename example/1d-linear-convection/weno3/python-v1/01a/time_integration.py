# time_integration.py

"""
时间推进器模块 (已集成注册系统)
使用装饰器 @register_component 自动注册
"""

from abc import ABC, abstractmethod
from cfd_registry import ComponentRegistry, register_component
from residual import ResidualCalculator
from boundary import BoundaryConditionFactory
            
# ---------------------- 1. 抽象时间推进器基类（统一接口） ----------------------
class TimeIntegrator(ABC):
    """时间推进器抽象基类：定义一维CFD时间推进的核心接口"""
    def __init__(self, cfd):
        self.cfd = cfd  # 持有CFD实例，获取配置/域/求解数据
        self.config = cfd.config
        self.domain = cfd.domain
        self.solution = cfd.solution
        self.residual_calculator = ResidualCalculator(cfd)
        self.boundary_condition = BoundaryConditionFactory.create(cfd)

    @abstractmethod
    def step(self, dt):
        """
        单次时间步推进（核心接口）
        :param dt: 时间步长
        :return: None
        """
        pass

    def init_total_field(self):
        self.solution.initialize_from_config(self.config)
        self.boundary_condition.apply(self.solution.u)
        self.solution.update_old_field()
    
    def compute_residual(self):
        """计算残差（所有RK方法都需要，封装为公共方法）"""
        self.residual_calculator.compute()

    def apply_boundary(self):
        """应用边界条件（公共逻辑）"""
        self.boundary_condition.apply(self.solution.u)

    def map_idx(self, i):
        """物理网格索引 → 残差数组索引（公共映射逻辑）"""
        return i - self.domain.ist

# ---------------------- 2. 具体RK时间推进器实现（使用装饰器注册） ----------------------

@register_component('integrator', 'rk1')
class RK1Integrator(TimeIntegrator):
    """1阶显式欧拉（RK1）"""
    def step(self, dt):
        self.compute_residual()
        for i in range(self.domain.ist, self.domain.ied):
            j = self.map_idx(i)
            self.solution.u[i] += dt * self.solution.res[j]
        self.apply_boundary()
        self.solution.update_old_field()
        
@register_component('integrator', 'rk2')
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
        
@register_component('integrator', 'rk3')
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

class TimeIntegratorFactory:
    """时间推进器工厂"""
    
    @staticmethod
    def create(cfd) -> 'TimeIntegrator':
        """
        创建时间推进器实例
        
        Args:
            cfd: CFD对象
        
        Returns:
            时间推进器实例
        """
        from factories.base_factory import BaseFactory
        rk_order = cfd.config.rk_order
        integrator_name = f'rk{rk_order}'
        return BaseFactory.create_component('integrator', integrator_name, cfd)