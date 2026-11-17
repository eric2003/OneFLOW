# core/registry.py
"""
统一注册系统 + 通用工厂
- ComponentRegistry: 组件注册表
- register_component: 装饰器
- BaseFactory: 通用工厂类
"""

from typing import Dict, Type, Any


# ==================== 1. 注册表核心 ====================
class ComponentRegistry:
    """组件注册表"""
    _registries: Dict[str, Dict[str, Type]] = {}
    _verbose = True

    @classmethod
    def set_verbose(cls, verbose: bool):
        cls._verbose = verbose

    @classmethod
    def register(cls, category: str, name: str, component_class: Type):
        #print(f"ComponentRegistry register cls={cls}")
        #print(f"ComponentRegistry register name={name}")
        #print(f"ComponentRegistry register component_class={component_class}")
        if category not in cls._registries:
            cls._registries[category] = {}
        if name in cls._registries[category]:
            if cls._registries[category][name] != component_class and cls._verbose:
                print(f"⚠️  覆盖注册: {category}.{name}")
        cls._registries[category][name] = component_class
        if cls._verbose:
            print(f"✅ 已注册: {category}.{name} -> {component_class.__name__}")

    @classmethod
    def get(cls, category: str, name: str) -> Type:
        if category not in cls._registries:
            raise ValueError(f"❌ 未知类别: {category} (可用: {list(cls._registries.keys())})")
        if name not in cls._registries[category]:
            raise ValueError(f"❌ 未找到: {category}.{name} (可用: {list(cls._registries[category].keys())})")
        return cls._registries[category][name]

    @classmethod
    def create(cls, category: str, name: str, *args, **kwargs) -> Any:
        component_class = cls.get(category, name)
        return component_class(*args, **kwargs)

    @classmethod
    def list_all(cls):
        return {cat: list(comps.keys()) for cat, comps in cls._registries.items()}


# ==================== 2. 装饰器 ====================
def register_component(category: str, name: str = None):
    """简化注册的装饰器"""
    def decorator(cls):
        component_name = name or cls.__name__.lower()
        #print(f"register_component  decorator name={name}")
        ComponentRegistry.register(category, component_name, cls)
        return cls
    return decorator


# ==================== 3. 通用工厂 ====================
class BaseFactory:
    """通用工厂基类"""
    @classmethod
    def create_component(cls, category: str, name: str, *args, **kwargs) -> Any:
        name_lower = name.lower()
        try:
            return ComponentRegistry.create(category, name_lower, *args, **kwargs)
        except ValueError as e:
            available = ComponentRegistry.list_all().get(category, [])
            if available:
                error_msg = f"不支持的 {category} 类型 '{name}'。可用类型：{available}"
            else:
                error_msg = f"不支持的 {category} 类型 '{name}'（无已注册组件）"
            raise ValueError(error_msg) from e

    @classmethod
    def get_available_components(cls, category: str) -> list:
        return ComponentRegistry.list_all().get(category, [])