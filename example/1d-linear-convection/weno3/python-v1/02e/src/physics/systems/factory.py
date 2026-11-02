# src/physics/systems/factory.py
"""
物理系统工厂
"""

from core.registry import BaseFactory

class PhysicalSystemFactory:
    """物理系统工厂"""
    
    @staticmethod
    def create(system_type: str, **kwargs):
        """
        创建物理系统实例
        
        Args:
            system_type: 系统类型，如 'linear_advection'
            **kwargs: 传递给方程的参数
        
        Returns:
            物理系统实例
        """
        return BaseFactory.create_component('system', system_type, **kwargs)
    
    @staticmethod
    def get_available_types():
        """获取所有可用的系统类型"""
        return BaseFactory.get_available_components('system')