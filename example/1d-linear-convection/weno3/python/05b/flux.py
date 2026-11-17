"""
通量计算器模块 (已集成注册系统)
使用装饰器 @register_component 自动注册，替代硬编码工厂
"""

from abc import ABC, abstractmethod

# ---------------------- 导入注册系统 ----------------------
try:
    # 从核心位置或本地导入注册系统
    from cfd_core.registry import ComponentRegistry, register_component
except ImportError:
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from registry import ComponentRegistry, register_component
        print("[flux] 使用本地 registry.py 中的注册系统")
    except ImportError:
        # 极简回退方案（实际使用时应该确保registry.py存在）
        print("[flux] 警告: 未找到注册系统，使用极简回退方案")
        class ComponentRegistry:
            _registries = {}
            @classmethod
            def register(cls, *args, **kwargs): pass
            @classmethod
            def create(cls, category, name, *args, **kwargs):
                # 回退到硬编码逻辑（仅用于过渡）
                if category == 'flux' and name == 'rusanov':
                    return RusanovFluxCalculator(*args, **kwargs)
                elif category == 'flux' and name == 'engquist-osher':
                    return EngquistOsherFluxCalculator(*args, **kwargs)
                raise ValueError(f"不支持的 {category}.{name}")
        def register_component(*args, **kwargs):
            def dummy_decorator(cls):
                return cls
            return dummy_decorator

# ---------------------- 1. 抽象通量计算基类（统一接口） ----------------------
class InviscidFluxCalculator(ABC):
    """无粘通量计算抽象基类：定义一维CFD通量计算接口"""
    def __init__(self, cfd):
        self.cfd = cfd
        self.config = cfd.config
        self.mesh = cfd.domain.mesh
        self.wave_speed = self.config.wave_speed

    @abstractmethod
    def compute(self, q_face_left, q_face_right, flux):
        """
        计算无粘通量（核心接口）
        :param q_face_left: 左界面值数组
        :param q_face_right: 右界面值数组
        :param flux: 输出通量数组
        :return: None
        """
        pass

# ---------------------- 2. 具体通量计算子类（使用装饰器注册） ----------------------

@register_component('flux', 'rusanov')
class RusanovFluxCalculator(InviscidFluxCalculator):
    """Rusanov（Lax-Friedrichs）通量"""
    def compute(self, q_face_left, q_face_right, flux):
        for i in range(self.mesh.nnodes):
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            c_L = self.wave_speed
            c_R = self.wave_speed
            F_L = c_L * u_L  # Flux from left state
            F_R = c_R * u_R  # Flux from right state
            Smax = max(abs(c_L), abs(c_R))  # Maximum wave speed
            flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L)
            
@register_component('flux', 'engquist-osher')
class EngquistOsherFluxCalculator(InviscidFluxCalculator):
    """Engquist-Osher通量（线性对流专用）"""
    def compute(self, q_face_left, q_face_right, flux):
        for i in range(self.mesh.nnodes):
            c = self.wave_speed
            cp = 0.5 * (c + abs(c))
            cm = 0.5 * (c - abs(c))
            u_L = q_face_left[i]
            u_R = q_face_right[i]
            flux[i] = cp * u_L + cm * u_R

# ---------------------- 3. 通量计算器工厂（使用注册系统） ----------------------
class FluxCalculatorFactory:
    """通量计算器工厂：根据配置创建对应通量计算器实例"""
    @staticmethod
    def create(cfd):
        """根据配置创建通量计算器实例"""
        flux_type = cfd.config.flux_type.lower()
        
        try:
            # 使用注册系统创建实例
            flux_instance = ComponentRegistry.create('flux', flux_type, cfd)
            return flux_instance
        except (ValueError, RuntimeError) as e:
            # 增强错误信息
            available = []
            try:
                all_comps = ComponentRegistry.list_all()
                available = all_comps.get('flux', [])
            except:
                # 后备列表
                available = ['rusanov', 'engquist-osher']
            
            raise ValueError(
                f"不支持的flux类型：'{flux_type}'\n"
                f"可用类型：{available}\n"
                f"原始错误：{e}"
            )