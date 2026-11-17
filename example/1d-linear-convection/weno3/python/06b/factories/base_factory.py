# factories/base_factory.py
"""
工厂基类 - 提供统一的组件创建接口
"""

from cfd_registry import ComponentRegistry
from typing import Any, Optional

class BaseFactory:
    """工厂基类，提供统一的组件创建方法"""
    
    @classmethod
    def create_component(cls, category: str, name: str, *args, **kwargs) -> Any:
        """
        统一创建组件的方法
        
        Args:
            category: 组件类别（如'boundary', 'flux'等）
            name: 组件名称（如'periodic', 'rusanov'等）
            *args, **kwargs: 传递给构造函数的参数
        
        Returns:
            创建的组件实例
        
        Raises:
            ValueError: 如果组件不存在
        """
        name_lower = name.lower()
        
        try:
            return ComponentRegistry.create(category, name_lower, *args, **kwargs)
        except ValueError as e:
            # 获取可用组件列表
            available = ComponentRegistry.list_all().get(category, [])
            
            # 构建友好的错误信息
            if available:
                error_msg = (f"不支持的{category}类型：'{name}'\n"
                           f"可用类型：{available}")
            else:
                error_msg = (f"不支持的{category}类型：'{name}'\n"
                           f"（{category}类别下没有注册任何组件）")
            
            raise ValueError(error_msg) from e
    
    @classmethod
    def get_available_components(cls, category: str) -> list:
        """获取指定类别的可用组件列表"""
        return ComponentRegistry.list_all().get(category, [])