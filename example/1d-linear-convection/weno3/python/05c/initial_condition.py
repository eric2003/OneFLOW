# initial_condition.py

"""
初始条件模块 (已集成注册系统)
使用装饰器 @register_component 自动注册
"""

import numpy as np
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
        print("[initial_condition] 使用本地 registry.py 中的注册系统")
    except ImportError:
        print("[initial_condition] 警告: 未找到注册系统，使用极简回退方案")
        class ComponentRegistry:
            _registries = {}
            @classmethod
            def register(cls, *args, **kwargs): pass
            @classmethod
            def create(cls, category, name, *args, **kwargs):
                # 回退到硬编码逻辑
                if category == 'initial_condition':
                    if name == 'step':
                        return StepFunctionIC(*args, **kwargs)
                    elif name == 'sin':
                        return SineWaveIC(*args, **kwargs)
                    elif name == 'gaussian':
                        return GaussianPulseIC(*args, **kwargs)
                raise ValueError(f"不支持的 {category}.{name}")
        def register_component(*args, **kwargs):
            def dummy_decorator(cls):
                return cls
            return dummy_decorator
            

# ---------------------- 1. 初始条件抽象基类 ----------------------
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


# ---------------------- 2. 具体初始条件实现（使用装饰器注册） ----------------------

@register_component('initial_condition', 'step')
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

@register_component('initial_condition', 'sin')
class SineWaveIC(InitialCondition):
    def evaluate_at(self, x):
        L = self.config.get("domain_length", 2.0)
        return np.sin(2 * np.pi * x / L)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)

@register_component('initial_condition', 'gaussian')
class GaussianPulseIC(InitialCondition):
    def evaluate_at(self, x):
        center = self.config.get("pulse_center", 0.5)
        width = self.config.get("pulse_width", 0.1)
        return np.exp(-((x - center) / width) ** 2)

    def apply(self, solution):
        x = solution.domain.mesh.xcc
        u0 = self.evaluate_at(x)
        self._apply_to_interior(solution, u0)


# ---------------------- 3. 初始条件工厂（使用注册系统） ----------------------
class InitialConditionFactory:
    """初始条件工厂：根据配置创建对应初始条件实例"""
    
    @classmethod
    def create(cls, ic_type, config):
        """创建初始条件实例"""
        ic_type_lower = ic_type.lower()
        
        try:
            # 使用注册系统创建实例
            ic_instance = ComponentRegistry.create('initial_condition', ic_type_lower, config)
            return ic_instance
        except (ValueError, RuntimeError) as e:
            # 增强错误信息
            available = []
            try:
                all_comps = ComponentRegistry.list_all()
                available = all_comps.get('initial_condition', [])
            except:
                # 后备列表
                available = ['step', 'sin', 'gaussian']
            
            raise ValueError(
                f"未知的初始条件类型: '{ic_type}'\n"
                f"支持的类型: {available}\n"
                f"原始错误: {e}"
            )
    
    # 保持原有的注册方法，用于向后兼容或动态注册
    _registry = {}  # 旧式注册表，可能被其他代码使用
    
    @classmethod
    def register(cls, name, ic_class):
        """注册初始条件类（兼容旧接口）"""
        cls._registry[name] = ic_class
        # 同时注册到新系统（如果可能）
        try:
            ComponentRegistry.register('initial_condition', name, ic_class)
        except:
            pass