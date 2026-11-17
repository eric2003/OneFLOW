# src/initial_conditions/factory.py
from core.registry import BaseFactory
from .base import InitialCondition

class InitialConditionFactory:
    """初始条件工厂"""
    
    @classmethod
    def create(cls, config) -> 'InitialCondition':
        """
        创建初始条件实例
        
        Args:
            config: 配置对象，包含ic_type属性
        
        Returns:
            初始条件实例
        """
        return BaseFactory.create_component('initial_condition', config.ic_type, config)
    
    @classmethod
    def get_available_types(cls):
        """获取所有可用的初始条件类型"""
        return BaseFactory.get_available_components('initial_condition')