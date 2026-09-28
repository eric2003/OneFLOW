# cfd_registry.py
"""
CFD注册系统统一导入模块
所有其他模块从这里导入注册系统
"""

import sys
import os

# 添加当前目录到路径，确保能找到registry.py
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    # 尝试导入注册系统
    from registry import ComponentRegistry, register_component
    REGISTRY_AVAILABLE = True
except ImportError:
    # 如果找不到，提供简单的空实现
    print("⚠️  警告: 未找到 registry.py，使用最小化回退实现")
    
    class ComponentRegistry:
        """最小化的注册系统回退实现"""
        _registries = {}
        
        @classmethod
        def register(cls, category, name, component_class):
            if category not in cls._registries:
                cls._registries[category] = {}
            cls._registries[category][name] = component_class
        
        @classmethod
        def get(cls, category, name):
            if category not in cls._registries:
                raise ValueError(f"未知类别: {category}")
            if name not in cls._registries[category]:
                raise ValueError(f"未找到: {category}.{name}")
            return cls._registries[category][name]
        
        @classmethod
        def create(cls, category, name, *args, **kwargs):
            component_class = cls.get(category, name)
            return component_class(*args, **kwargs)
        
        @classmethod
        def list_all(cls):
            return {cat: list(comps.keys()) 
                   for cat, comps in cls._registries.items()}
    
    def register_component(category, name=None):
        """简化的装饰器"""
        def decorator(cls):
            component_name = name or cls.__name__.lower()
            ComponentRegistry.register(category, component_name, cls)
            return cls
        return decorator
    
    REGISTRY_AVAILABLE = False

# 导出统一的接口
__all__ = ['ComponentRegistry', 'register_component', 'REGISTRY_AVAILABLE']