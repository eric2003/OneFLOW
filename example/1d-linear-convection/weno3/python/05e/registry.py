"""
CFD组件注册系统核心
完全独立，不依赖任何现有代码
"""

from typing import Dict, Type, Any

class ComponentRegistry:
    """组件注册表 - 替代硬编码工厂"""
    
    # 存储所有注册的组件 {类别: {名称: 类}}
    _registries: Dict[str, Dict[str, Type]] = {}
    
    @classmethod
    def register(cls, category: str, name: str, component_class: Type):
        """注册一个组件类"""
        if category not in cls._registries:
            cls._registries[category] = {}
        
        cls._registries[category][name] = component_class
        print(f"✅ 已注册: {category}.{name} -> {component_class.__name__}")
    
    @classmethod
    def get(cls, category: str, name: str) -> Type:
        """获取已注册的组件类"""
        if category not in cls._registries:
            raise ValueError(f"❌ 未知类别: {category}")
        
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


def register_component(category: str, name: str = None):
    """简化注册的装饰器"""
    def decorator(cls):
        component_name = name or cls.__name__.lower().replace('boundary', '').replace('flux', '')
        ComponentRegistry.register(category, component_name, cls)
        return cls
    return decorator