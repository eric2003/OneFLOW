# registry.py (改进版)
"""
CFD组件注册系统 - 简洁高效版
"""

from typing import Dict, Type, Any

class ComponentRegistry:
    """组件注册表"""
    
    _registries: Dict[str, Dict[str, Type]] = {}
    _verbose = True  # 控制是否打印注册信息
    
    @classmethod
    def set_verbose(cls, verbose: bool):
        """设置是否打印注册信息"""
        cls._verbose = verbose
    
    @classmethod
    def register(cls, category: str, name: str, component_class: Type):
        """注册组件类"""
        if category not in cls._registries:
            cls._registries[category] = {}
        
        # 检查是否已注册
        if name in cls._registries[category]:
            if cls._registries[category][name] == component_class:
                # 相同类重复注册，静默跳过
                return
            elif cls._verbose:
                print(f"⚠️  覆盖注册: {category}.{name}")
        
        cls._registries[category][name] = component_class
        if cls._verbose:
            print(f"✅ 已注册: {category}.{name} -> {component_class.__name__}")
    
    @classmethod
    def get(cls, category: str, name: str) -> Type:
        """获取组件类"""
        if category not in cls._registries:
            available = list(cls._registries.keys())
            raise ValueError(f"❌ 未知类别: {category} (可用类别: {available})")
        
        if name not in cls._registries[category]:
            available = list(cls._registries[category].keys())
            raise ValueError(f"❌ 未找到: {category}.{name} (可用: {available})")
        
        return cls._registries[category][name]
    
    @classmethod
    def create(cls, category: str, name: str, *args, **kwargs) -> Any:
        """创建组件实例"""
        component_class = cls.get(category, name)
        return component_class(*args, **kwargs)
    
    @classmethod
    def list_all(cls) -> Dict[str, list]:
        """列出所有已注册的组件"""
        return {
            category: list(components.keys())
            for category, components in cls._registries.items()
        }
    
    @classmethod
    def get_count(cls) -> int:
        """获取注册组件总数"""
        return sum(len(components) for components in cls._registries.values())

def register_component(category: str, name: str = None):
    """简化注册的装饰器"""
    def decorator(cls):
        component_name = name or cls.__name__.lower()
        ComponentRegistry.register(category, component_name, cls)
        return cls
    return decorator