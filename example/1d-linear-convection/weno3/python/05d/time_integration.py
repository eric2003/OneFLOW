# time_integration.py

"""
时间推进器模块 (已集成注册系统)
使用装饰器 @register_component 自动注册
"""

from abc import ABC, abstractmethod

# ---------------------- 导入注册系统 ----------------------
try:
    from cfd_core.registry import ComponentRegistry, register_component
except ImportError:
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from registry import ComponentRegistry, register_component
        print("[time_integration] 使用本地 registry.py 中的注册系统")
    except ImportError:
        print("[time_integration] 警告: 未找到注册系统，使用极简回退方案")
        class ComponentRegistry:
            _registries = {}
            @classmethod
            def register(cls, *args, **kwargs): pass
            @classmethod
            def create(cls, category, name, *args, **kwargs):
                # 回退到硬编码逻辑
                if category == 'integrator':
                    if name == 'rk1':
                        return RK1Integrator(*args, **kwargs)
                    elif name == 'rk2':
                        return RK2Integrator(*args, **kwargs)
                    elif name == 'rk3':
                        return RK3Integrator(*args, **kwargs)
                raise ValueError(f"不支持的 {category}.{name}")
        def register_component(*args, **kwargs):
            def dummy_decorator(cls):
                return cls
            return dummy_decorator
            
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

# ---------------------- 3. 时间推进器工厂（使用注册系统） ----------------------
class TimeIntegratorFactory:
    """时间推进器工厂：根据配置创建对应RK实例"""
    
    @staticmethod
    def create(cfd):
        """根据配置创建时间推进器实例"""
        rk_order = cfd.config.rk_order
        
        # 将数字顺序映射为注册系统使用的名字
        # 注意：注册时使用了 'rk1', 'rk2', 'rk3' 这样的名字
        integrator_name = f'rk{rk_order}'
        
        try:
            # 使用注册系统创建实例
            integrator_instance = ComponentRegistry.create('integrator', integrator_name, cfd)
            return integrator_instance
        except (ValueError, RuntimeError) as e:
            # 增强错误信息
            available = []
            try:
                all_comps = ComponentRegistry.list_all()
                available = all_comps.get('integrator', [])
            except:
                # 后备列表
                available = ['rk1', 'rk2', 'rk3']
            
            # 尝试提供更友好的错误信息
            available_orders = [int(name[2:]) for name in available if name.startswith('rk')]
            
            raise ValueError(
                f"不支持的RK阶数：{rk_order}\n"
                f"支持的RK阶数：{available_orders}\n"
                f"原始错误：{e}"
            )